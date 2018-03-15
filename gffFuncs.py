import pprint
from BCBio.GFF import GFFExaminer

in_file = "mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20161111.gff"
examiner = GFFExaminer()
in_handle = open(in_file)
pprint.pprint(examiner.parent_child_map(in_handle))
in_handle.close()
