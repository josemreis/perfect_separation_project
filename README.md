# Separation project

Term project for a class on *Analysis of Discrete data* at Nova IMS.

1. Clone the repo

```bash
git clone https://github.com/josemreis/separation_project.git
```

2. Open R in the main directory of this folder and render the poster

   ```R
   ### packages
   packs <- c("rmarkdown", "posterdown")
   for (pack in packs) {
     
     if (pack %in% installed.packages()){
       
       require(pack, character.only = TRUE)
       
     } else {
       
       install.packages(pack, dependencies = TRUE)
       require(pack, character.only = TRUE)
       
     }  
     
   }
   ### render
   rmarkdown::render("./separation_poster/separation_poster.Rmd")
   ```

   