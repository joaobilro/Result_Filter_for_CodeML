# !/usr/bin/env python3
#
# -*- coding: utf-8 -*-
#
# Site_Models_CMLFR.py
#
# This Python 3 script allows the user to filter the results in the output files created 
# with CodeML from the Phylogenetic Analysis by Maximum Likelihood (PAML) package by Z. Yang (1997), 
# namely by running the site models analyses, on a set of user-defined genes.
# In order for it to work as intended, arrange the Null and Alternative runs for each gene as follows:
#
# /Main_Directory/
#         /Gene_A/
#             /Null/
#                <Contents of the null CodeML site models run for gene A>
#             /Alternative/
#                 <Contents of the alternative CodeML site models run for gene A>
#         /Gene_B/
#             /Null/
#                <Contents of the null CodeML site models run for gene B>
#             /Alternative/
#                 <Contents of the alternative CodeML site models run for gene B>
#         /Gene_C/
#             /Null/
#                 <Contents of the null CodeML site models run for gene C>
#             /Alternative/
#                 <Contents of the alternative CodeML site models run for gene C>
#        ...
#
# It is imperative that for each gene's folder, the Null and Alternative subfolders be present as previously
# exemplified, otherwise the script will not work. The name for each gene's folder does not matter as long as
# it does not conflict "basic naming rules" (no special characters, blank spaces, etc.).
# Also, the output files should have an unique extension, such as .mlc or .otp. This can be altered in the script.
#           _
#         ><_> 
#     
#        
# MIT License
# 
# Copyright (c) 2024 João Bilro (joaobilro)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import os

cmlfr = argparse.ArgumentParser(description="This Python 3 script allows the user to filter the results in the output files created"
                                 "with CodeML from the Phylogenetic Analysis by Maximum Likelihood (PAML) package by Z. Yang (1997)," 
                                 "namely by running the site models analyses, on a set of user-defined genes.")

cmlfr.add_argument("--input", "-i", dest="input_gene_list", required=True, nargs="+", help="The main directory containing all the structured gene folders.")

cmlfr.add_argument("--output", "-o", dest="output_file", required=True, help="The output .csv file, containing the filtered results.")

args = cmlfr.parse_args()

class SiteModelsResults:
    """ This class contains the output values for both Null and Alternative runs (and their comparisons) for the user-specified genes
    contained in a given main directory. If the LRT is significant, further information for each gene will be available in the ouput."""

    def 