# PeacoQC flow cytometry QC operator

#### Description

`peacoqc` operator performs quality control on flow cytometry data.

##### Usage

Input projection|.
---|---
`row`   | represents the variables (e.g. channels, markers)
`col`   | represents the observations (Use 'Time' on top of rowid.) 
`y-axis`| measurement value


Output relations|.
---|---
`QC_flag`       | character, quality flag `pass` or `fail`

Input parameters|.
---|---
`MAD`         | The MAD parameter. Default is 6. If this is increased, the algorithm becomes less strict
`IT_limit`         | The IsolationTree parameter. Default is 0.55. If this is increased, the algorithm becomes less strict
`remove_zeros`     | Default is FALSE. If TRUE, the zero values will be removed before the peak detection step. They will not be indicated as ’bad’ value. Recommended when cleaning mass cytometry data

##### Details
The peacoQC operator will determine peaks on the channels in the flow cytometry data. 
Then it will flagg anomalies caused by e.g. clogs, changes in speed etc. by using an IsolationTree and the MAD method (median absolute deviation). The parameters can be changed to make the quality checks more or less strict.
The operator returns a quality flagg `pass` or `fail`. 


#### Reference

[PeacoQC R package]((http://www.bioconductor.org/packages/release/bioc/html/PeacoQC.html))

##### See Also

[flowAI operator](https://github.com/tercen/flowai_operator)
[flowClean operator](https://github.com/tercen/flowclean_operator)
[flowCut operator](https://github.com/tercen/flowcut_operator)
