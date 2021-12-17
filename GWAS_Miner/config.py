# NLP Variables
regex_entity_patterns = {
    "PVAL": [
        r"(?:\(?[pP][  =<-]+(?:val[ue]*)?)(\d*\.?.?\d*[ ]?[*×xX][ ]?\d+([ (]?[-−]\d+[ ]?)*[-−]\d+)",
        # r"(?:[pP][  =<-]){1,}(?:val){0,}(?:ue){0,}([  <≥=×xX-]{0,}[  \(]?\d+[\.]?[\d]{0,}[-−_^*()  \d×xX]{0,})"
        r"(?:[pP][  =<-])+(?:val)*(?:ue)*(?:[  =-])*([<≥=×xX]*(?:[  \(])*\d+[\.]?[\d]*(?:[(])?[-−_^*  \d×xX]*)"
        # r"(?:[pP][VAL]{0,3})[ =<>]+([0-9]{1,}[\.]?[0-9]{0,})",
        # r"(\d?\..?\d[ ]?[*×xX]{1}[ ]?\d{1,}[ (]?[-−]\d{1,}[ )]?)",
        # r"((\(?\b[pP][  =<-]{1,}(val{1,}[ue]{0,})?[  <≥=×xX-]{0,}[  \(]?\d+[\.]?[\d]{0,}[-^*()  \d×xX]{0,"
        # r"})|(\d?\.?\d[  ]?[*×xX]{1}[  ]?\d{1,}[  ]?-\d{1,}))",
        # "([pP][- ]{1,2}[val]{0,3}[ue]{0,}[^0-9]{0,9})([0-9]+)([0-9.e-]+)",
        # r"([pP][VAL]{0,3}[ =]+[xX× _\-−]+[0-9]+)"
    ],
    # "PTYPE": r"(\(?GEE\)?)|(\(?FBAT\)?)",
    # "Table Ref": r"(table[- ]{0,}\d{1,})"
}

pheno_assoc_patterns = [
    [
        {
            "RIGHT_ID": "phenotype_anchor",
            "RIGHT_ATTRS": {
                "_": {
                    "is_trait": True
                }
            }
        },
        {
            "LEFT_ID": "phenotype_anchor",
            "REL_OP": ".*",
            "RIGHT_ID": "marker_subject",
            "RIGHT_ATTRS": {"DEP": {"IN": ["dep", "nsubj","dobj", "conj", "amod", "compound", "appos", "parataxis"]}, "ENT_TYPE": "RSID"},
        },
        {
            "LEFT_ID": "marker_subject",
            "REL_OP": "$++",
            "RIGHT_ID": "significance_subject",
            "RIGHT_ATTRS": {"DEP": {"IN": ["dep", "nsubj","dobj", "conj", "amod", "compound", "appos"]}, "ENT_TYPE": "PVAL"},
        }
    ],
    [
        {
            "RIGHT_ID": "phenotype_anchor",
            "RIGHT_ATTRS": {
                "_": {
                    "is_trait": True
                }
            }
        },
        {
            "LEFT_ID": "phenotype_anchor",
            "REL_OP": ".*",
            "RIGHT_ID": "marker_subject",
            "RIGHT_ATTRS": {"DEP": {"IN": ["dep", "nsubj","dobj", "conj", "amod", "compound", "appos", "parataxis"]}, "ENT_TYPE": "RSID"},
        },
        {
            "LEFT_ID": "marker_subject",
            "REL_OP": "$--",
            "RIGHT_ID": "significance_subject",
            "RIGHT_ATTRS": {"DEP": {"IN": ["dep", "nsubj","dobj", "conj", "amod", "compound", "appos"]}, "ENT_TYPE": "PVAL"},
        }
    ],
    [
        {
            "RIGHT_ID": "phenotype_anchor",
            "RIGHT_ATTRS": {
                "_": {
                    "is_trait": True
                }
            }
        },
        {
            "LEFT_ID": "phenotype_anchor",
            "REL_OP": ";*",
            "RIGHT_ID": "marker_subject",
            "RIGHT_ATTRS": {"DEP": {"IN": ["dep", "nsubj","dobj", "conj", "amod", "compound", "appos", "parataxis"]}, "ENT_TYPE": "RSID"},
        },
        {
            "LEFT_ID": "marker_subject",
            "REL_OP": "$++",
            "RIGHT_ID": "significance_subject",
            "RIGHT_ATTRS": {"DEP": {"IN": ["dep", "nsubj","dobj", "conj", "amod", "compound", "appos"]}, "ENT_TYPE": "PVAL"},
        }
    ],
    [
        {
            "RIGHT_ID": "phenotype_anchor",
            "RIGHT_ATTRS": {
                "_": {
                    "is_trait": True
                }
            }
        },
        {
            "LEFT_ID": "phenotype_anchor",
            "REL_OP": ";*",
            "RIGHT_ID": "marker_subject",
            "RIGHT_ATTRS": {"DEP": {"IN": ["dep", "nsubj","dobj", "conj", "amod", "compound", "appos", "parataxis"]}, "ENT_TYPE": "RSID"},
        },
        {
            "LEFT_ID": "marker_subject",
            "REL_OP": "$--",
            "RIGHT_ID": "significance_subject",
            "RIGHT_ATTRS": {"DEP": {"IN": ["dep", "nsubj", "dobj", "conj", "amod", "compound", "appos"]}, "ENT_TYPE": "PVAL"},
        }
    ]
]
