%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE REQUIREMENTS %%%%
%
% [1] GEM/Data/common_tasks_growth_RPMI1640.xlsx
%     Common growth task list.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%
%
% [1] iBrain2845/Results/iBrain2845.sbml.xlsx
%     Generalised GEM.
%
% [2] main/Results/DESeq2/
%     DESeq2 results.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUTS %%%%
%
% [1] main/Results/GEM/repMets_cX_vs_cY.tsv
%     Reporter metabolite analysis results between clusters X and Y.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%

% Import
refmodel = importExcelModel('iBrain2845/Results/iBrain2845.sbml.xlsx', false);
refmodel.rxnNames = refmodel.rxns;
nClust = 4;
phenotypes(1) = "Alzheimer's disease";
phenotypes(2) = "Control";
phenotypes(3) = "MCI";
phenotypes(4) = "Other dementia";
infile_DESeq2 = 'main/Results/DESeq2/';

% Reporter metabolite analysis
%
for i = 1:nClust
    for j = 1:nClust
        if i ~= j
            repMet.infile = strcat(infile_DESeq2, 'allDEG_', convertStringsToChars(phenotypes(i)), '_vs_', convertStringsToChars(phenotypes(j)), '.tsv');
            repMet.DESeq2 = readTXT(repMet.infile);
            repMet.genes = repMet.DESeq2(2:end,1);
            repMet.l2FC = str2double(repMet.DESeq2(2:end,3));
            repMet.pvalues = str2double(repMet.DESeq2(2:end,6));
            repMet.outfile = strcat(outfile_dir, 'repMets_', convertStringsToChars(phenotypes(i)), '_vs_', convertStringsToChars(phenotypes(j)), '.txt');
            repMet.repMets = reporterMetabolites(refmodel, repMet.genes, repMet.pvalues, true, repMet.outfile, repMet.l2FC);
        end
    end
end
