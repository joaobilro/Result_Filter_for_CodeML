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
from collections import OrderedDict
import os
from scipy.stats import chi2
import sys

cmlfr = argparse.ArgumentParser(description="This Python 3 script allows the user to filter the results in the output files created"
                                 "with CodeML from the Phylogenetic Analysis by Maximum Likelihood (PAML) package by Z. Yang (1997)," 
                                 "namely by running the site models analyses, on a set of user-defined genes.")

cmlfr.add_argument("--input", "-i", dest="input_gene_list", required=True, nargs="+", help="The main directory containing all the structured gene folders.")

cmlfr.add_argument("--output", "-o", dest="output_file", required=True, help="The output .csv file, containing the filtered results.")

args = cmlfr.parse_args()

class SiteModelsResults ():
    """Contains the output values for both Null and Alternative runs (and their comparisons) for the user-specified genes present in
    a given main directory. If the LRT is significant, further information for each gene will be available in the ouput."""

    def extract_null_values(main_directory):
        """Extract all the values and information needed from the output files of CodeML's null run."""
    
    ### This will list all the folders that start with "gene_", which is the adopted structure and naming
    gene_folders = [dir for dir in os.listdir() if os.path.isdir(os.path.join(os.getcwd(), dir)) and dir.startswith("gene_")]
    
    for folder in gene_folders:
        ### This will establish the path to the Null folder inside each gene folder
        output_path = os.path.join(os.getcwd(), folder, "Null")
        for file in os.listdir(output_path):
            ### Get the path to the first (and only) .mlc file in the Null folder
            if file.endswith(".mlc"):
                output_file_path = os.path.join(output_path, output_file_path)
                print("Null: Found .mlc file in", folder)
                break
            else:
                print("Error: No .mlc file found in", output_path)
                sys.exit(1) ### Abort the script if a .mlc file is missing     
        
        ### Now, to open the file and start parsing the information
        output_file = open(output_file_path, 'r')
        
        for line in output_file:
            ### Check if the first line is empty, if not it contains the number of seqs/individuals followed by the gene length 
            if line.strip().startswith("Before deleting alignment gaps"):
                ### Extract gene length from the next line after the term
                gene_length = next(output_file).strip().split()[1]
                
            elif line.strip().startswith("After deleting gaps."):
                ### Extract gene length from the next line after the term
                gene_length = next(output_file).strip().split()[1]
            
            ### If none of the previous terms are present in the file...
            elif line.strip():
                ### Extract gene length from the first line
                gene_length = int(line.split()[1])
    
            else:
                print("Error: Gene length not found in", output_file_path)
                sys.exit(1) ### Abort the script if there is a .mlc file without the gene length
            
            if line.strip().startswith("lnL"):
                m0_lnL = float(line.split(":")[-1].split()[0])
            else:
                print("Error: Null lnL not found in", output_file_path)
                sys.exit(1) ### Abort the script if there is a .mlc file without the lnL
        
            if line.strip().startswith("omega (dN/dS)"):
                omega = float(line.split("=")[1].strip())
            else:
                print("Error: Omega not found in", output_file_path)
                sys.exit(1) ### Abort the script if there is a .mlc file without omega (dN/dS)
                            
        output_file.close()
        
    
    def extract_alternative_values(main_directory):
        """Extract all the values and information needed from the output files of CodeML's alternative run."""
    
    ### This will list all the folders that start with "gene_", which is the adopted structure and naming
    gene_folders = [dir for dir in os.listdir() if os.path.isdir(os.path.join(os.getcwd(), dir)) and dir.startswith("gene_")]
    
    for folder in gene_folders:
        ### This will establish the path to the Alternative folder inside each gene folder
        output_path = os.path.join(os.getcwd(), folder, "Alternative")
        for file in os.listdir(output_path):
            ### Get the path to the first (and only) .mlc file in the Alternative folder
            if file.endswith(".mlc"):
                output_file_path = os.path.join(output_path, output_file_path)
                print("Alternative: Found .mlc file in", folder)
                break
            else:
                print("Error: No .mlc file found in", output_path)
                sys.exit(1) ### Abort the script if a .mlc file is missing     

        ### Now, to open the file and start parsing the information
        output_file = open(output_file_path, 'r')
        
        for line in output_file:
            ### Check if the first line is empty, if not it contains the number of seqs/individuals followed by the gene length 
            if line.strip().startswith("Before deleting alignment gaps"):
                ### Extract gene length from the next line after the term
                gene_length = next(output_file).strip().split()[1]
                
            elif line.strip().startswith("After deleting gaps."):
                ### Extract gene length from the next line after the term
                gene_length = next(output_file).strip().split()[1]
            
            ### If none of the previous terms are present in the file...
            elif line.strip():
                ### Extract gene length from the first line
                gene_length = int(line.split()[1])
    
            else:
                print("Error: Gene length not found in", output_file_path)
                sys.exit(1) ### Abort the script if there is a .mlc file without the gene length
            
            ### Counter to find the correct lnL values
            lnl_counter = 0
            
            if line.strip().startswith("lnL"):
                ### First lnL is for the null model
                lnl_counter += 1
                
                if lnl_counter == 2: 
                    ### Second occurrence is for the M1a model
                    m1a_lnl = float(line.split(":")[-1].split()[0])
                
                if lnl_counter == 3:
                    ### Third occurrence is for the M2a model
                    m2a_lnl = float(line.split(":")[-1].split()[0])
                    
                    if line.strip().startswith("Bayes Empirical Bayes (BEB) analysis"):
                        
                
                if lnl_counter == 4:
                    ### Fourth occurrence is for the M7 model
                    m7_lnl = float(line.split(":")[-1].split()[0])
                
                if lnl_counter == 5:
                    ### Fifth occurrence is for the M8 model
                    m8_lnl = float(line.split(":")[-1].split()[0]) 
            else:
                print("Error: Null lnL not found in", output_file_path)
                sys.exit(1) ### Abort the script if there is a .mlc file without the lnL
        
        output_file.close()                
                    
    def likelihood_ratio_test(df, null_lnl, alt_lnl):
        """Perform a likelihood ratio test (LRT) with the lnL values for the Null and Alternative models, and then calculate the p-value using a Chi-squared test."""
        
        ### Need both lnL values to be able to perform a LRT
        if alt_lnl is None or null_lnl is None:
            lrt = 0
            pvalue = 1.0
        else:
            lrt = 2 * (alt_lnl - (null_lnl))
            ### Negative LRTs do not exist, it is the same as having LRT = 0, which is the lower limit for LRT
            if lrt <= 0:
                pvalue = 1.0
            else:
                ### Calculating pvalue with scipy.stats
                pvalue = float(1 - chi2(df).cdf(lrt))        
        
        return lrt, pvalue
                
        
if __name__ == "__main__":
    def main():
        ### Matching arguments with their intended variables
        main_directory = args.input_gene_list
        output_csv = args.output_file
            
        