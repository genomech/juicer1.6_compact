import pysam
import json
import bisect
import subprocess

def LoadFragmentMap(RestrSitesMap):
	FragmentMap = {}
	with open(RestrSitesMap, 'rt') as MapFile:
		for Contig in MapFile:
			List = Contig[:-1].split(' ')
			FragmentMap[List[0]] = [int(item) for item in List[1:]]
	return FragmentMap

def CalcDist(Item1, Item2): 
	if ((Item1 is None) or (Item2 is None) or (type(Item1) == list) or (type(Item2) == list)): return None
	return float("+inf") if (Item1["Chr"] != Item2["Chr"]) else abs(Item1["Pos"] - Item2["Pos"])

def SortItems(Item1, Item2): return tuple([(item["ID"], item["Pos"]) for item in sorted([Item1, Item2], key=lambda x: (x["RefID"], x["Pos"]))])

def ProcessQuery(Query, ChromSizes, MinMAPQ):
	# Filter unmapped
	if any([item[1].is_unmapped for item in Query["ReadBlock"]]): return { "ReadBlock": Query["ReadBlock"], "Type": "Unmapped" }
	if any([item[1].mapping_quality < MinMAPQ for item in Query["ReadBlock"]]): return { "ReadBlock": Query["ReadBlock"], "Type": "MappingQualityFailed" }
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
		Type = ("1" if item.is_read1 else "2") + ("s" if (item.is_secondary or item.is_supplementary) else "p")
		TypeDict[Type].append({ "ID": int(index), "Chr": str(item.reference_name), "RefID": int(item.reference_id), "Pos": int(End) if item.is_reverse else int(Start) })
	# Create Pattern
	Pattern = tuple([len(item) for index, item in TypeDict.items()])
	TypeDict = { index: (None if not item else (item[0] if len(item) == 1 else item)) for index, item in TypeDict.items() }
	Dist = { f"1{index1}2{index2}": CalcDist(TypeDict[f"1{index1}"], TypeDict[f"2{index2}"]) for index1, index2 in ('pp', 'ps', 'sp', 'ss')}
	# Norm Chimera 4 Ends
	if Pattern == (1, 1, 1, 1):
		if ((Dist["1p2p"] < 1000) and (Dist["1s2s"] < 1000)) or ((Dist["1p2s"] < 1000) and (Dist["1s2p"] < 1000)):
			Sorted = SortItems(TypeDict["1p"], TypeDict["1s"])
			Pair = [{ "Read": Query["ReadBlock"][Sorted[0][0]][1], "Pos": Sorted[0][1] }, { "Read": Query["ReadBlock"][Sorted[1][0]][1], "Pos": Sorted[1][1] }]
			return { "ReadBlock": Query["ReadBlock"], "Type": "ChimericPaired", "Pair": Pair }
		else: return { "ReadBlock": Query["ReadBlock"], "Type": "ChimericAmbiguous" }
	# Norm Chimera 3 Ends
	elif Pattern in ((1, 0, 1, 1), (1, 1, 1, 0)):
		if TypeDict["1s"] is None:
			if ((Dist["1p2p"] < 1000) or (Dist["1p2s"] < 1000)): Sorted = SortItems(TypeDict["1p"], TypeDict["2p"] if Dist["1p2p"] > Dist["1p2s"] else TypeDict["2s"])
			else: Sorted = None
		if TypeDict["2s"] is None:
			if ((Dist["1p2p"] < 1000) or (Dist["1s2p"] < 1000)): Sorted = SortItems(TypeDict["2p"], TypeDict["1p"] if Dist["1p2p"] > Dist["1s2p"] else TypeDict["1s"])
			else: Sorted = None
		if Sorted is None: return { "ReadBlock": Query["ReadBlock"], "Type": "ChimericAmbiguous" }
		Pair = [{ "Read": Query["ReadBlock"][Sorted[0][0]][1], "Pos": Sorted[0][1] }, { "Read": Query["ReadBlock"][Sorted[1][0]][1], "Pos": Sorted[1][1] }]
		return { "ReadBlock": Query["ReadBlock"], "Type": "ChimericPaired", "Pair": Pair }
	# Regular Pair
	elif Pattern == (1, 0, 1, 0):
		Sorted = SortItems(TypeDict["1p"], TypeDict["2p"])
		Pair = [{ "Read": Query["ReadBlock"][Sorted[0][0]][1], "Pos": Sorted[0][1] }, { "Read": Query["ReadBlock"][Sorted[1][0]][1], "Pos": Sorted[1][1] }]
		return { "ReadBlock": Query["ReadBlock"], "Type": "NormalPaired", "Pair": Pair }
	# Collisions
	elif (Pattern[1] > 1) or (Pattern[3] > 1):
		pass # TODO Collisions
	# Other
	return { "ReadBlock": Query["ReadBlock"], "Type": "ChimericAmbiguous" }

def Main(InputFileSAM, OutputFileTXT, InterPairsTXT, ChimericAmbiguousFileSAM, UnmappedSAM, MappingQualityFailedSAM, StatsTXT, RestrictionSiteFile = None, MinMAPQ = 0):
	Input = pysam.AlignmentFile(InputFileSAM, 'r', check_sq=False)
	SortCommand = (f'sort -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n | tee >( gzip -c > "{OutputFileTXT}" ) |' +
				f' awk -F " " \'{{print $2 "\\t" $3 "\\t" $6 "\\t" $7}}\' | gzip -c > "{InterPairsTXT}"')
	Output = subprocess.Popen(SortCommand, shell=True, executable="/bin/bash", stdin=subprocess.PIPE)
	if RestrictionSiteFile is not None: FragmentMap = LoadFragmentMap(RestrictionSiteFile)
	TechInfo = {
		"ChimericAmbiguous": pysam.AlignmentFile(ChimericAmbiguousFileSAM, "wb", template = Input),
		"Unmapped": pysam.AlignmentFile(UnmappedSAM, "wb", template = Input),
		"MappingQualityFailed": pysam.AlignmentFile(UnmappedSAM, "wb", template = Input)
		}
	ChromSizes = { Input.references[i]: Input.lengths[i] for i in range(Input.nreferences) }
	Stats = { "SequencedReadPairs": 0, "NormalPaired": 0, "ChimericPaired": 0, "ChimericAmbiguous": 0, "MappingQualityFailed": 0, "Unmapped": 0, "Ligation": { "Motif": None, "LineCount": 0, "PresentCount": 0 } }
	Query = { "ReadName": None, "ReadBlock": [] }
	
	def BlockProcess():
		Stats["SequencedReadPairs"] += 1
		Query["ReadBlock"] = list(enumerate(Query["ReadBlock"]))
		Result = ProcessQuery(Query, ChromSizes, MinMAPQ)
		Stats[Result["Type"]] += 1
		if Result["Type"] in ("Unmapped", "ChimericAmbiguous", "MappingQualityFailed"):
			for index, Rec in Query["ReadBlock"]: TechInfo[Result["Type"]].write(Rec)
		if Result["Type"] in ("ChimericPaired", "NormalPaired"):
			Read1, Read2 = Result["Pair"]
			Line = ' '.join([
				'16' if Read1["Read"].is_reverse else '0',
				str(Read1["Read"].reference_name),
				str(Read1["Pos"]),
				'0' if RestrictionSiteFile is None else str(bisect.bisect(FragmentMap[Read1["Read"].reference_name], Read1["Pos"])),
				'16' if Read2["Read"].is_reverse else '0',
				str(Read2["Read"].reference_name),
				str(Read2["Pos"]),
				'1' if RestrictionSiteFile is None else str(bisect.bisect(FragmentMap[Read2["Read"].reference_name], Read2["Pos"])),
				str(Read1["Read"].mapping_quality),
				str(Read1["Read"].cigarstring),
				str(Read1["Read"].seq.__str__()),
				str(Read2["Read"].mapping_quality),
				str(Read2["Read"].cigarstring),
				str(Read2["Read"].seq.__str__()),
				str(Read1["Read"].query_name),
				str(Read2["Read"].query_name)
				]) + '\n'
			Output.stdin.write(Line.encode('utf-8'))
		
	while 1:
		try:
			Record = next(Input)
			if not (Record.is_secondary or Record.is_supplementary):
				Stats["Ligation"]["LineCount"] += 1
				# TODO Add ligation counter
			if Record.query_name == Query["ReadName"]: Query["ReadBlock"].append(Record)
			else:
				BlockProcess()
				Query["ReadName"] = Record.query_name
				Query["ReadBlock"].clear()
				Query["ReadBlock"].append(Record)
		except StopIteration:
			BlockProcess()
			Input.close()
			Output.stdin.close()
			Output.wait()
			Stats["Alignable"] = Stats["ChimericPaired"] + Stats["NormalPaired"]
			for stat in ("ChimericPaired", "ChimericAmbiguous", "NormalPaired", "Unmapped", "Alignable", "MappingQualityFailed"): Stats[stat] = { "Count": Stats[stat], "%": Stats[stat] / Stats["SequencedReadPairs"] * 100 }
			Stats["Ligation"]["%"] = Stats["Ligation"]["PresentCount"] / Stats["SequencedReadPairs"] * 100 # BUG WTF?
			# TODO Postprocessing? Library Complexity?
			json.dump(Stats, open(StatsTXT, 'wt'), indent=4, ensure_ascii=False)
			break

Main(InputFileSAM = "/Data/NGS_Data/20211228_NGS_MinjaF_Pool/Results/Human_HiC/K1/splits/8_S73_L003.fastq.gz.filtered.sam", OutputFileTXT = "test_mergednodups.txt.gz", InterPairsTXT = "test_interpairs.txt.gz", MappingQualityFailedSAM = "/dev/null", ChimericAmbiguousFileSAM = "/dev/null", UnmappedSAM = "/dev/null", StatsTXT = "test.stats.txt", RestrictionSiteFile = None, MinMAPQ = 30)
