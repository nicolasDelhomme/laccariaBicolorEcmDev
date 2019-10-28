# laccariaBicolorEcmDev
The molecular timecourse of hybrid aspen/Laccaria bicolor ectomycorrhiza development

## Abstract

The project aims at deciphering the molecular timecourse of the establishment of Laccaria bicolor ectomycorrhiza on hybrid aspen (Populus tremula x Populus tremuloides, T89). Samples collected from various stage of the ectomycorrhiza have been collected and sequenced

## Setup
### Repository
In the terminal, clone both project and its submodules
```{bash git,eval=FALSE}
git clone --recurse-submodules git@github.com:nicolasDelhomme/laccariaBicolorEcmDev.git  
```

Or in two steps:
```{bash git submodule,eval=FALSE}
git clone git@github.com:nicolasDelhomme/laccariaBicolorEcmDev.git
cd laccariaBicolorEcmDev
git submodule init
git submodule update
```
### Data
```{bash setup,eval=FALSE}
ln -s /mnt/picea/projects/aspseq/jfelten/T89-Laccaria-bicolor data
```

### Update
 To update the submodule when doing `git pull`, you need to run extra commands:
 
 ```{bash git update,eval=FALSE}
 cd laccariaBicolorEcmDev
 git submodule update --remote
 ```
 
 Note that this will overwrite anyt local changes, also commited, by the content of the [UPSCb-common repository](git@github.com:UPSCb/UPSCb-common.git)
 
 Therefore any changes to be done to these common scripts is best done directly in the 
 [UPSCb-common repository](git@github.com:UPSCb/UPSCb-common.git)
 
 