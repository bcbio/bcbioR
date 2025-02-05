Using the brown fat scRNA-seq dataset (TN vs cold7), we ran both propeller and sccomp to evaluate the results. For a majority of the celltypes there as concordance showing that proprtions of cells were significantly different between conditions. Howver, for celltypes identified by only one of the methods the result was not always convincing. Best practice for the template may be to run both methods and compare results. The most conservative approach would be to take the intersect of the two methods.

<p align="center">
<img src="comp.png" width="500">
</p>


Another method, currently being used in the training material, is MiloR. This tool transforms to neighborhoods so is not directly comparable to propeller and sccomp.



