"""An OOP library! OMG! Wraps various other executables to make I/O easier (or
at least I.)

External programs wrapped:
HYPHYMD (MacOS), version ???
"""

class HyphyRelax:
    """Wrapper class to call HYPHY RELAX. So far, only supports 1 input file
    containing both the sequence and tree data; doesn't support those two bits
    as separate files.
    """
    def __init__(self):
        self.analysis_type="1" # Selection analysis
        self.analysis_subtype="7" # RELAX
        self.genetic_code="1" # Universal genetic code
        self.data_file="" # input data file, with both sequence and tree data.
        self.use_tree="y" # Use the tree file in the input data_file.
        self.branch_test_set="2" # Use labelled branches. 1 = Use unlabelled.
        self.relax_analysis_type="2" # Run only 2 model. 1 = run 4 models.

    def get_cmd(self, outfile_suffix="_results.out"):
        '''Get a single bash command string. There's probably a better way to
        do this.
        '''
        cmd1 = "(echo "+self.analysis_type+"; "
        cmd2 = "echo "+self.analysis_subtype+"; "
        cmd3 = "echo "+self.genetic_code+"; "
        cmd4 = "echo "+self.data_file+"; "
        cmd5 = "echo "+self.use_tree+"; "
        cmd6 = "echo "+self.branch_test_set+"; "
        cmd7 = "echo "+self.relax_analysis_type+"; "
        cmd_suffix = "| HYPHYMP > "+self.data_file+outfile_suffix

        cmd = cmd1 + cmd2 + cmd3 + cmd4 + cmd5 + cmd6 + cmd7 + cmd_suffix

        return cmd

    def show_me(self):
        """Human Helper: Displays all current attributes and values of the
        HyphyRelax object.
        """
        print("analysis_type = %s" % self.analysis_type)
        print("analysis_subtype = %s" % self.analysis_subtype)
        print("genetic_code = %s" % self.genetic_code)
        print("data_file = %s" % self.data_file)
        print("use_tree = %s" % self.use_tree)
        print("branch_test_set = %s" % self.branch_test_set)
        print("relax_analysis_type = %s" % self.relax_analysis_type)
