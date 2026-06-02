run "full_data_analysis" in the parent dataset with the estrela dataset. The section that produces the coupling matrix should 
output all panels in SI Fig.22.

Beware: by default, the code takes out the relative-abundance mode and produces panel B, and all subpanels in C except for PC14(=LC1).
To plot te interaction network including the relative abundance mode, set the parameter "startindex=0" in the coupling matrix section.

Figures are outputted to "out/Estrella_full/Figures/"

The individual eigenvectors can be produced from the "plot_eigenvector" script and will be outputted to "individual_eigenvector"
