from marshmallow import fields, Schema, validates_schema, ValidationError


class Param(Schema):
    chr_col = fields.Int(required=True, description="Column number for chromosome")
    pos_col = fields.Int(required=True, description="Column number for base position")
    ea_col = fields.Int(required=True, description="Column number for effect allele")
    oa_col = fields.Int(required=True, description="Column number for other allele")
    beta_col = fields.Int(required=True, description="Column number for effect")
    se_col = fields.Int(required=True, description="Column number for standard error")
    pval_col = fields.Int(required=True, description="Column number for association P value")
    delimiter = fields.Str(required=True, description="Input file column delimiter")
    header = fields.Bool(required=True, description="Does the input file have a header")
    ncase_col = fields.Int(required=False, description="Column number for number of cases")
    snp_col = fields.Int(required=False, description="Column number for variant identifier")
    eaf_col = fields.Int(required=False,
                         description="Column number for effect allele variant frequency")
    oaf_col = fields.Int(required=False,
                         description="Column number for other allele variant frequency")
    imp_z_col = fields.Int(required=False,
                           description="Column number for summary statistics imputation Z score")
    imp_info_col = fields.Int(required=False,
                              description="Column number for summary statistics imputation INFO score")
    ncontrol_col = fields.Int(required=False,
                              description="Column number for number of controls (if case/control) or total sample size if continuous")
    build = fields.Str(required=True, description="Name of the genome build i.e. GRCh36, GRCh37, GRCh38")
    data = fields.Str(required=False,
                      description="Path to input file, required if not passed in main command line parameters")
    out = fields.Str(required=False,
                     description="Path to output vcf file, required if not passed in main command line parameters")
    id = fields.Str(required=False,
                    description="Unique identified or name of GWAS study, required if not passed in main command line parameters")
    cohort_cases = fields.Float(required=False,
                                description="Total study number of cases")
    cohort_controls = fields.Float(required=False,
                                   description="Total study number of controls (if case/control) or total sample size if continuous")
    jsonmeta = fields.Str(required=False, description="Path to input metadata file file")
    md5 = fields.Str(required=False, description="md5 checksum for input file")

    @validates_schema(pass_original=True)
    def check_unknown_fields(self, data, original_data):
        unknown = set(original_data) - set(self.fields)
        if unknown:
            raise ValidationError('Unknown field', unknown)
