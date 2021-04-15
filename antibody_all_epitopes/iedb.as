table iedb
"iedb bed9+ track"
   (
   string chrom;       "Reference sequence chromosome or scaffold"
   uint   chromStart;  "Start position in chromosome"
   string name;        "Name or ID of item, ideally both human readable and unique"
   uint score;         "Score (0-1000)"
   char[1] strand;     "+ or - for strand"
   uint reserved;       	"RGB value (use R,G,B string in input file)"
   string epitopeID;    "Epitope ID in IEDB"
   string epitopeType;  "Linear or Non-linear epitope"
   string aaStartPos;     "Amino acid start position"
   string aaEndPos;       "Amino acid end position"
   string antigenName;  "Antigen name"
   string antigenAccession; "Antigen accession"
   string organismName; "Organism Name"
   string assayType; "Assay Type"
   string antibodyName; "Antibody Name"
   )
