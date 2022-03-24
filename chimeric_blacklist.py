import pysam
import json

def CalcDist(Item1, Item2): 
	if ((Item1 is None) or (Item2 is None)): return None
	return float("+inf") if (Item1["Chr"] != Item2["Chr"]) else abs(Item1["Pos"] - Item2["Pos"])

def SortItems(Item1, Item2): return tuple([(item["ID"], item["Pos"]) for item in sorted([Item1, Item2], key=lambda x: (x["RefID"], x["Pos"]))])

def ProcessQuery(Query, ChromSizes):
	# Filter unmapped
	if any([item[1].is_unmapped for item in Query["ReadBlock"]]): return { "ReadBlock": Query["ReadBlock"], "Type": "Unmapped" }
	# Create Sorter
	TypeDict = { index: list() for index in ("1p", "1s", "2p", "2s") }
	# Annotation
	for index, item in Query["ReadBlock"]:
		Start = item.reference_start + 1
		End = item.reference_end
		CigarFirst = item.cigar[0]
		CigarLast = item.cigar[-1]
		SoftHard = (4, 5)
		if CigarFirst[0] in SoftHard:
			Start -= CigarFirst[1]
			if Start <= 0: Start = 1
		if CigarLast[0] in SoftHard:
			End += CigarLast[1]
			if End >= ChromSizes[item.reference_name]: End = ChromSizes[item.reference_name]
		Type = ("1" if item.is_read1 else "2") + ("s" if item.is_secondary else "p")
		TypeDict[Type].append({ "ID": int(index), "Chr": str(item.reference_name), "RefID": int(item.reference_id), "Pos": int(End) if item.is_reverse else int(Start) })
	# Create Pattern
	Pattern = tuple([len(item) for index, item in TypeDict.items()])
	TypeDict = { index: None if not item else item[0] for index, item in TypeDict.items() }
	Dist = { f"1{index1}2{index2}": CalcDist(TypeDict[f"1{index1}"], TypeDict[f"2{index2}"]) for index1, index2 in ('pp', 'ps', 'sp', 'ss')}
	# Norm Chimera 4 Ends
	if Pattern == (1, 1, 1, 1):
		if ((Dist["1p2p"] < 1000) and (Dist["1s2s"] < 1000)) or ((Dist["1p2s"] < 1000) and (Dist["1s2p"] < 1000)):
			Sorted = SortItems(TypeDict["1p"], TypeDict["1s"])
			Pair = [{ "Read": Query["ReadBlock"][Sorted[0][0]][1], "Pos": Sorted[0][1] }, { "Read": Query["ReadBlock"][Sorted[1][0]][1], "Pos": Sorted[1][1] }]
			return { "ReadBlock": Query["ReadBlock"], "Type": "NormChimera", "Pair": Pair }
		else: return { "ReadBlock": Query["ReadBlock"], "Type": "AbnormChimera" }
	# Norm Chimera 3 Ends
	elif Pattern in ((1, 0, 1, 1), (1, 1, 1, 0)):
		if TypeDict["1s"] is None:
			if ((Dist["1p2p"] < 1000) or (Dist["1p2s"] < 1000)): Sorted = SortItems(TypeDict["1p"], TypeDict["2p"] if Dist["1p2p"] > Dist["1p2s"] else TypeDict["2s"])
			else: Sorted = None
		if TypeDict["2s"] is None:
			if ((Dist["1p2p"] < 1000) or (Dist["1s2p"] < 1000)): Sorted = SortItems(TypeDict["2p"], TypeDict["1p"] if Dist["1p2p"] > Dist["1s2p"] else TypeDict["1s"])
			else: Sorted = None
		if Sorted is None: return { "ReadBlock": Query["ReadBlock"], "Type": "AbnormChimera" }
		Pair = [{ "Read": Query["ReadBlock"][Sorted[0][0]][1], "Pos": Sorted[0][1] }, { "Read": Query["ReadBlock"][Sorted[1][0]][1], "Pos": Sorted[1][1] }]
		return { "ReadBlock": Query["ReadBlock"], "Type": "NormChimera", "Pair": Pair }
	# Regular Pair
	elif Pattern == (1, 0, 1, 0):
		Sorted = SortItems(TypeDict["1p"], TypeDict["2p"])
		Pair = [{ "Read": Query["ReadBlock"][Sorted[0][0]][1], "Pos": Sorted[0][1] }, { "Read": Query["ReadBlock"][Sorted[1][0]][1], "Pos": Sorted[1][1] }]
		return { "ReadBlock": Query["ReadBlock"], "Type": "RegularPairs", "Pair": Pair }
	# Other
	return { "ReadBlock": Query["ReadBlock"], "Type": "AbnormChimera" }

def Main(InputFileSAM, OutputFileTXT, AbnormChimeraFileSAM, UnmappedSAM, StatsTXT):
	Input = pysam.AlignmentFile(InputFileSAM, 'r', check_sq=False)
	Output = open(OutputFileTXT, 'wt')
	TechInfo = {
		"AbnormChimera": pysam.AlignmentFile(AbnormChimeraFileSAM, "wb", template = Input),
		"Unmapped": pysam.AlignmentFile(UnmappedSAM, "wb", template = Input)
		}
	ChromSizes = { Input.references[i]: Input.lengths[i] for i in range(Input.nreferences) }
	Stats = { "Total": 0, "NormChimera": 0, "AbnormChimera": 0, "RegularPairs": 0, "Unmapped": 0 }
	Query = { "ReadName": None, "ReadBlock": [] }
	
	def BlockProcess():
		Stats["Total"] += 1
		Query["ReadBlock"] = list(enumerate(Query["ReadBlock"]))
		Result = ProcessQuery(Query, ChromSizes)
		Stats[Result["Type"]] += 1
		if Result["Type"] in ("Unmapped", "AbnormChimera"):
			for index, Rec in Query["ReadBlock"]: TechInfo[Result["Type"]].write(Rec)
		if Result["Type"] in ("NormChimera", "RegularPairs"):
			Read1, Read2 = Result["Pair"]
			Line = [
				'16' if Read1["Read"].is_reverse else '0',
				str(Read1["Read"].reference_name),
				str(Read1["Pos"]),
				'16' if Read2["Read"].is_reverse else '0',
				str(Read2["Read"].reference_name),
				str(Read2["Pos"]),
				str(Read1["Read"].mapping_quality),
				str(Read1["Read"].cigarstring),
				str(Read1["Read"].seq.__str__()),
				str(Read2["Read"].mapping_quality),
				str(Read2["Read"].cigarstring),
				str(Read2["Read"].seq.__str__()),
				str(Read1["Read"].query_name),
				str(Read2["Read"].query_name)
				]
			Output.write(' '.join(Line) + '\n')
		
	while 1:
		try:
			Record = next(Input)
			if Record.query_name == Query["ReadName"]: Query["ReadBlock"].append(Record)
			else:
				BlockProcess()
				Query["ReadName"] = Record.query_name
				Query["ReadBlock"].clear()
				Query["ReadBlock"].append(Record)
		except StopIteration:
			BlockProcess()
			Input.close()
			Output.close()
			json.dump(Stats, open(StatsTXT, 'wt'), indent=4, ensure_ascii=False)
			break
