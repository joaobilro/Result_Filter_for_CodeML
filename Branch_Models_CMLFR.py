# !/usr/bin/env python3
#
# -*- coding: utf-8 -*-
#
# Branch_Models_CMLFR.py
#
# Version 1.0
#
# This Python 3 script allows the user to filter the results in the output files created 
# with CodeML from the Phylogenetic Analysis by Maximum Likelihood (PAML) package by Z. Yang (1997), 
# namely by running the branch models analyses, on a set of user-defined genes.
# In order for it to work as intended, arrange the Null and Alternative runs for each gene as follows:
#
# /Main_Directory/
#         /Gene_A/
#             /Null/
#                <Contents of the null CodeML branch models run for gene A>
#             /Alternative/
#                 <Contents of the alternative CodeML branch models run for gene A>
#         /Gene_B/
#             /Null/
#                <Contents of the null CodeML branch models run for gene B>
#             /Alternative/
#                 <Contents of the alternative CodeML branch models run for gene B>
#         /Gene_C/
#             /Null/
#                 <Contents of the null CodeML branch models run for gene C>
#             /Alternative/
#                 <Contents of the alternative CodeML branch models run for gene C>
#        ...
#
# It is imperative that for each gene's folder, the Null and Alternative subfolders be present as previously
# exemplified, otherwise the script will not work. The name for each gene's folder does not matter as long as
# it does not conflict "basic naming rules" (no special characters, blank spaces, etc.).
# This script does not allow to parse through different group's results at the same time. It should be individually
# executed for each group, inside each group's folder. Each of these folders should maintain the specified structure. 
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
import glob
import os
from scipy.stats import chi2
import sys

cmlfr = argparse.ArgumentParser(description="This Python 3 script allows the user to filter the results in the output files created"
                                 "with CodeML from the Phylogenetic Analysis by Maximum Likelihood (PAML) package by Z. Yang (1997)," 
                                 "namely by running the branch models analyses, on a set of user-defined genes.")

cmlfr.add_argument("--input", "-i", dest="input_gene_list", required=True, type = str, help="The main directory containing all the structured gene folders.")

cmlfr.add_argument("--output", "-o", dest="output_file", required=True, type = str, help="The output .csv file, containing the filtered results.")

args = cmlfr.parse_args()

class GeneFolders:
    """Contains the data and structure for each gene folder, which will be used to populate the .csv output file."""

    def __init__(self):
        
        self.folder_name = None                          ### The name for the gene folder
        self.gene_length = None                          ### The length of the gene
        self.omega = None                                ### The omega (dN/dS) estimated by CodeML's null run for the gene
        self.m0_lnl = None                               ### The log likelihood value for the M0 model
        self.m2a_lnl = None                              ### The log likelihood value for the M2a model
        self.df0v2a = 1                                  ### The degrees of freedom for the comparison between the M0 and M2a models
        self.lrt0v1a = None                              ### The LRT from the comparison between the M0 and M2a models
        self.pvalue0v2a = None                           ### The pvalue calculated by a Chi Squared test with the LRT value from the comparison between the M0 and M2a models
        self.foreground_omega = None                     ### The foreground omega (dN/dS) estimated by CodeML's alternative run for the gene
        self.background_omega = None                     ### The background omega (dN/dS) estimated by CodeML's alternative run for the gene

class BranchModelsResults:
    """Contains the output values for both Null and Alternative runs (and their comparisons) for the user-specified genes present in
    a given main directory. If the LRT is significant, further information for each gene will be available in the ouput."""

    def __init__(self, main_directory, output_csv):
        self.dict = OrderedDict()
        self.main_directory = main_directory
        self.output_csv = output_csv

    def extract_null_values(self):
        """Extract all the values and information needed from the output files of CodeML's null run."""
    
        ### This will list all the folders that start with "gene_", which is the adopted structure and naming
        os.chdir(self.main_directory)
        gene_folders = glob.glob("gene_*")

        mlc_found = False   ### This flag will be used to check if there is at least one .mlc file inside each run's folder

        for folder in gene_folders:
            values = GeneFolders()
            values.folder_name = folder.split("gene_")[1]
            ### This will establish the path to the Null folder inside each gene folder
            output_path = os.path.join(os.getcwd(), folder, "Null")
            file_list = os.listdir(output_path)
            
            for file in file_list:
                ### Get the path to the first (and only) .mlc file in the Null folder
                if file.endswith(".mlc"):
                    mlc_found = True
                    gene_length = -1
                    output_file_path = os.path.join(output_path, file)
                    print("Null: Found .mlc file in", folder)    
            
                    ### Now, to open the file and start parsing the information
                    output_file = open(output_file_path, "r")
                    
                    for line in output_file:
                        ### Check if the first line is empty, if not it contains the number of seqs/individuals followed by the gene length 
                        if gene_length == -1:
                            ### Extract gene length from the first line
                            split = line.split(" ")
                            gene_length = split[-1].split("\n")[0]
                            values.gene_length = gene_length
                        
                        ### If gaps/ambiguous sites are not preserved...
                        if line.strip().startswith("After deleting gaps."):
                            ### Extract gene length from the next line after the term
                            gene_length = next(output_file).strip().split()[1]
                            values.gene_length = gene_length
                            
                        if line.strip().startswith("lnL"):
                            m0_lnL = float(line.split(":")[-1].split()[0])
                            values.m0_lnl = m0_lnL
                        
                        if line.strip().startswith("omega (dN/dS)"):
                            omega = float(line.split("=")[1].strip())
                            values.omega = omega

                    output_file.close()

                    self.dict[values.folder_name] = values      
                
                    break   
            
            if not mlc_found:
                print("Error: No .mlc file found in", folder)
                sys.exit(1) ### Abort the script if there is a missing .mlc file

    def extract_alternative_values(self):
        """Extract all the values and information needed from the output files of CodeML's alternative run."""
    
        ### This will list all the folders that start with "gene_", which is the adopted structure and naming
        os.chdir(self.main_directory)
        gene_folders = glob.glob("gene_*")
        
        mlc_found = False   ### This flag will be used to check if there is at least one .mlc file inside each run's folder

        for folder in gene_folders:
            values = self.dict[folder.split("gene_")[1]]
            ### This will establish the path to the Alternative folder inside each gene folder
            output_path = os.path.join(os.getcwd(), folder, "Alternative")
            file_list = os.listdir(output_path)
            
            for file in file_list:
                ### Get the path to the first (and only) .mlc file in the Alternative folder
                if file.endswith(".mlc"):
                    mlc_found = True
                    output_file_path = os.path.join(output_path, file)
                    print("Alternative: Found .mlc file in", folder)    

                    ### Now, to open the file and start parsing the information
                    output_file = open(output_file_path, 'r')
                    

                    for line in output_file:
                        
                       if line.strip().startswith("lnL"):
                            m2a_lnL = float(line.split(":")[-1].split()[0])
                            values.m2a_lnl = m2a_lnL
                       
                       if line.strip().startswith("w (dN/dS) for branches"):
                           background = float(line.split(":")[-1].split()[0])
                           values.background_omega = background
                           foreground = float(line.split(":")[-1].split()[1])
                           values.foreground_omega = foreground


                    output_file.close()                

                    self.dict.update([(values.folder_name, values)])
                
                    break
            
            if not mlc_found:
                print("Error: No .mlc file found in", folder)
                sys.exit(1) ### Abort the script if there is a missing .mlc file

    def likelihood_ratio_test(self, df, null_lnl, alt_lnl):
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
    
    def df1_lrt(self):
        ### Comparison between M0 and M2a
        for gene in self.dict.values():
                gene.lrt0v2a, gene.pvalue0v2a = self.likelihood_ratio_test(gene.df0v2a, gene.m0_lnl , gene.m2a_lnl)   ### df = 1 

    def write_output_csv(self):
        """Write the contents of the ordered dictionary to an output .csv file."""            

        save_to_output = open(self.output_csv, "w")
        ### write header to output .csv
        save_to_output.write("Gene;Length (bp);Omega (dN/dS);M0 lnL;M2a lnL;df;LRT (0 vs 2a);pvalue;Foreground dN/dS;Background dN/dS\n")
        
        for gene in self.dict.values():
            print("Processing results for {}...".format(gene.folder_name))

            if gene.pvalue0v2a <= 0.05:
                
                save_to_output.write("{};{};{};{};{};{};{};{};{};{}\n".format(
                                                gene.folder_name,                      
                                                gene.gene_length,                      
                                                gene.omega,                                  
                                                gene.m0_lnl,                                
                                                gene.m2a_lnl,                              
                                                gene.df0v2a,                                
                                                gene.lrt0v2a,                              
                                                gene.pvalue0v2a,                        
                                                gene.foreground_omega,
                                                gene.background_omega))
            
            else:
                save_to_output.write("{};{};{};{};{};{};{};{};{};{}\n".format(
                                                gene.folder_name,                      
                                                gene.gene_length,                      
                                                gene.omega,                                  
                                                gene.m0_lnl,                                
                                                gene.m2a_lnl,                              
                                                gene.df0v2a,                                
                                                gene.lrt0v2a,                              
                                                gene.pvalue0v2a,                        
                                                None,
                                                None))

        save_to_output.close()

def main():
        ### Matching arguments with their intended variables
        main_directory = args.input_gene_list
        output_csv = args.output_file

        sitemodels = BranchModelsResults(main_directory, output_csv)

        ### Start running the functions
        sitemodels.extract_null_values()
        sitemodels.extract_alternative_values()
        sitemodels.df1_lrt()

        ### Save everything to the output file
        sitemodels.write_output_csv()

        ### Finished
        print("The script has finished parsing the results.")

if __name__ == "__main__":

    main()    