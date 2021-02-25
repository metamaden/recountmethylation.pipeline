#!/usr/bin/env Rscript

# Author: Sean Maden
#
# Preprocess available sample metadata
#

# require(data.table); require(rjson)

#-------------------------
# md postprocess functions
#-------------------------
# disease terms (general, non-cancer)
#' Disease terms for disease variable
#' 
#' List containing search terms and string queries for regex pattern matching.
#' 
#' @return List of terms (names) and string queries (values).
#' @export
md_post_disease <- function(){
  dxl <- list("acute" = c("acute"), "syndrome" = c("syndrome"), 
              "disorder" = c("disorder"), "inflam"= c("inflammation"), 
              "cancer" = c("cancer"), "case" = c("case"), 
              "normal" = c("normal"), 
              "healthy" = c("healthy"), "replicate" = c("replicate"),
              "control" = c("control", "CONTROL", "ctl", "CTL", "ctrl", 
                            "CTRL"),
              "psychosis" = c("psychosis", "psychotic"), 
              "schizophrenia" = c("schizophrenia"),
              "arthritis" = c("osteoarthritis", "arthritis", "rheumatoid", 
                              "psoriatic", "fibromyalgia", "gout"),
              "rheumatoid_arthritis" = c("rheumatoid arthritis"), 
              "osteoarthritis" = c("osteoarthritis"), 
              "psoriatic_arthritis" = c("psoriatic arthritis"),
              "genetic_disorder" = c("fragile x", "cystic fibrosis", "duane", 
                                     "polycystic", "chrons", "hemophelia", 
                                     "haemophelia", "hemochromatosis", 
                                     "huntington's", "huntingtons", 
                                     "thalassemia", "tay sachs", "tay sach", 
                                     "parkinson's", "parkinsons", 
                                     "sickle cell", "marfan"), 
              "parkinsons" = c("parkinsons", "parkinson's"), 
              "sickle_cell" = c("sickle cell"), 
              "anemia" = c("anemia", "sickle cell"),
              "alzheimers" = c("alzheimer", "alzheimer"),
              "dementia" = c("dementia"), "lewy_body" = c("lewy bod"),
              "cystic_fibrosis" = c("cystic fibrosis"), 
              "scoliosis" = c("scoliosis"), 
              "obese" = c("obese"), 
              "irritable_bowel_disease" = c("IBD", "ibd", "irritable bowel"),
              "lesions" = c("lesions"), 
              "insulin_resistance" = c("insulin resist"),
              "autism" = c("autism", "autistic"), "patient" = c("patient"))
  return(dxl)
}

# cancer tissue terms/generic cancer terms
#' Cancer terms for tissue variable
#' 
#' List containing search terms and string queries for regex pattern matching.
#' 
#' @return List of terms (names) and string queries (values).
#' @export
md_post_cancer <- function(){
  cxl <- list("cancer" = c("tumor", "tumour", "metasta", "carcinoma", 
                           "sarcoma", "neoplas", "adenoma", "cancroid"),
              "tumor" = c("tumor", "tumour"),"metastasis" = c("metasta"),
              "carcinoma" = c("carcinoma"),"sarcoma" = c("sarcoma"), 
              "neoplasia" = c("neoplas"), "adenoma" = c("adenoma"), 
              "adenocarcinoma" = c("adenocarcinoma"), "lepidic" = c("lepidic"),
              "blastoma" = c("blastoma"), "benign" = c("benign"))
  return(cxl)
}

# cancer types by tissue/location
#' Cancer terms list for tissue variable
#' 
#' List containing search terms and string queries for regex pattern matching.
#' 
#' @return List of terms (names) and string queries (values).
#' @export
md_post_cancertype <- function(){
  cxsubl <- list("skin_cancer" = c("melanoma", "skin cancer"),
                 "brain_cancer" = c("glioblastoma", "astrocytoma", 
                                    "brain cancer", 
                                    "medulloblastoma"),
                 "medulloblastoma" = c("medulloblastoma"), 
                 "glioblastoma" = c("glioblastoma"),
                 "breast_cancer" = c("breast lobular carcinoma", 
                                     "breast ductal carcinoma", 
                                     "breast cancer", "triple negative"),
                 "colorectal_cancer" = c("colorectal adeno", "colon cancer", 
                                         "colorectal cancer", "rectal cancer"),
                 "stomach_cancer" = c("stomach adeno", "stomach cancer", 
                                      "gastric cancer", "gastric adeno"),
                 "esophageal_cancer" = c("esophageal carcinoma", 
                                         "esophageal adeno",
                                         "esophageal squamous cell carcinoma",
                                         "oesophageal carcinoma", 
                                         "oesophageal adeno", 
                                         "oesophageal squamous cell carcinoma",
                                         "esophageal cancer", 
                                         "oesophageal cancer", 
                                         " EAC$"),
                 "nerve_cell_cancer" = c("paraganglioma", "ependymoma", 
                                         "nerve cancer", 
                                         "nerve cell cancer", 
                                         "schwannoma"),
                 "paraganglioma" = c("paraganglioma"),
                 "ovarian_cancer" = c("ovarian serous carcinoma", 
                                      "ovarian cancer", 
                                      "endometrioid", 
                                      "ovarian epithelial cancer"),
                 "uterine_cancer" = c("uterine carcinosarcoma", 
                                      "uterine cancer",
                                      "uterine corpus endometrial carcinoma", 
                                      "endometrial carcinoma", 
                                      "uterine serous carcinoma", 
                                      "uterine papillary serous carcinoma"),
                 "kidney_cancer" = c("oncocytoma", 
                                     "clear cell renal cell carcinoma", 
                                     "chromophobe renal cell carcinoma", 
                                     "renal cancer", "kidney cancer",
                                     "kidney papillary carcinoma"),
                 "thyroid_cancer" = c("thyroid carcinoma", "thyroid cancer"),
                 "lung_cancer" = c("lung adenocarcinoma", 
                                   "lung squamous cell carcinoma", 
                                   "lung cancer", 
                                   "non-mucinous bronchoalveolar carcinoma",
                                   "lepidic-predominant adenocarcinoma", 
                                   "LPA"),
                 "bladder_cancer" = c("invasive urothelial bladder cancer", 
                                      "bladder cancer"),
                 "prostate_cancer" = c("prostate adenocarcinoma", 
                                       "prostate cancer"),
                 "liver_cancer" = c("liver hepatocellular carcinoma", 
                                    "hepatoblastoma", 
                                    "cholangiocarcinoma", 
                                    "liver angiosarcoma", 
                                    "liver cancer"),
                 "thymus_gland_cancer" = c("thymoma", "thymus cancer", 
                                           "thymus gland cancer", 
                                           "thymic cancer"),
                 "testicular_cancer" = c("testicular germ cell cancer", 
                                         "testicular cancer"),
                 "pancreatic_cancer" = c("pancreatic ductal adenocarcinoma", 
                                         "pancreatic cancer"),
                 "cervical_cancer" = c("cervical squamous cell carcinoma", 
                                       "cervical cancer",
                                       "cervical squamous cell adenocarcinoma"
                                       ),
                 "eye_cancer" = c("uveal melanoma", "uveal lymphoma", 
                                  "intraocular cancer", "retinoblastoma",
                                  "retinal cancer"))
  return(cxsubl)
}

# leukemia disease terms 
#' Leukemia terms list for disease variable
#' 
#' List containing search terms and string queries for regex pattern matching.
#' 
#' @return List of terms (names) and string queries (values).
#' @export
md_post_leukemia <- function(){
  leukl <- list("leukemia" = c("leukemia", "chronic leuk", "chronic myelo", 
                               "acute leuk", "acute lympho", "acute myel", 
                               "cml", "CML", "aml", "AML", "ALL"),
                "acute_leukemia" = c("acute leuk","acute lympho", "acute myel",
                                     "aml", "AML", "ALL"),
                "acute_myeloid_leukemia" = c("acute myel", "aml", "AML"),
                "acute_lymphoblastic_leukemia" = c("acute lympho", "ALL"))
  return(leukl)
}

# tissue terms, blood and related
#' Blood terms list for tissue variable
#' 
#' List containing search terms and string queries for regex pattern matching.
#' 
#' @return List of terms (names) and string queries (values).
#' @export
md_post_tissue_blood <- function(){
  txl <- list("blood" = c("blood", "hematopoiet", "haematopoiet","lymphoid", 
                          "myeloid", "natural killer", "( |^)nk( |$)", 
                          "( |^)NK( |$)", "erythrocyte", "mast cell", 
                          "myeloblast", "plasma","monocyte", "lymphocyte", 
                          "eosinophil", "neutrophil","basophil", "macrophage",
                          "megakaryocyte", "thrombocyte","wbc", "WBC", "rbc", 
                          "RBC","bcell", "b cell", "tcell", "t cell", "cd4", 
                          "cd5", "cd8", "cd34", "CD4", "CD5", "CD8", "CD34", 
                          "cytotoxic", "helper", "peripheral blood leukocytes", 
                          "( |^)PBL( |$|:)", "pbmc", "( |^)PBMC( |$|:)", 
                          "buffy", "blood spot", "blood punch", "granulocyte",
                          "white blood cell"),
              "buffy_coat" = c("buffy"), 
              "whole_blood" = c("whole blood"), 
              "peripheral_blood" = c("peripheral blood"), 
              "cord_blood" = c("cord blood"), "blood_spot" = c("blood spot"), 
              "white_blood_cell" = c("wbc", "WBC", "white blood cell", 
                                     "monocyte", "lymphocyte", "eosinophil", 
                                     "neutrophil", "basophil","bcell", 
                                     "b cell", "tcell", "t cell", "cd4", 
                                     "cd5", "cd8", "cd34", "granulocyte", 
                                     "CD4", "CD5", "CD8", "CD34", "cytotoxic", 
                                     "helper", "pbmc", "PBMC"),
              "peripheral_blood_mononuclear_cells" = c("pbmc", "PBMC"), 
              "peripheral_blood_leukocytes" = c("peripheral blood leukocytes", 
                                                "( |^)PBL( |$|:)"), 
              "cd4" = c("cd4", "CD4"), "cd5" = c("cd5", "CD5"), 
              "cd8" = c("cd8", "CD8"), "cd34" = c("cd34", "CD34"), 
              "granulocyte" = c("granulocyte"), "monocyte" = c("monocyte"), 
              "lymphocyte" = c("lymphocyte"), "neutrophil" = c("neutrophil"), 
              "eosinophil" = c("eosinophil"), "basophil" = c("basophil"), 
              "t_cell" = c("tcell", "t cell", "cd4", "cd5", "cd8", 
                           "cd34", "cytotoxic", "helper"))
  return(txl)
}

# tissue terms, general, not blood
#' Non-blood terms list for tissue variable
#' 
#' List containing search terms and string queries for regex pattern matching.
#' 
#' @return List of terms (names) and string queries (values).
#' @export
md_post_tissue <- function(){
  which.var <- c("gsm_title", "sample_type")
  txl <- list("adjacent" = c("adjacent"), "matched" = c("match"), 
              "distal" = c("distal"), "medullary" = c("medullary"), 
              "cultured" = c("cultured"), "paired" = c("paired"), 
              "explant" = c("explant"), "biopsy" = c("biopsy"), 
              "clone" = c("clone", "clonal"), 
              "subclone" = c("subclone", "subclonal"),
              "resection" = c("resection"), "xenograft" = c("xenograft"), 
              "cells" = c("cells"), "cell_line" = c("cell line"), 
              "peripheral" = c("peripheral"), "whole" = c("whole"),
              "organoid" = c("organoid"), 
              "parenchyma" = c("parenchyma", "parenchyme"), 
              "mesenchyme" = c("mesenchyme", "mesoderm"), 
              "mesoderm" = c("mesoderm"), 
              "ectoderm" = c("ectoderm"), "endoderm" = c("endoderm"),
              "gland" = c("gland"), "acinus" = c("acinus", "acinary"),
              "muscle" = c("muscle", "brachii", "ulnaris", "minimus",
                           "gemellus", "gluteus", "bicep", "rhomboids",
                           "tongue", "splenius", "capitis"),
              "smooth_muscle" = c("smooth muscle", "smoothe muscle"), 
              "skeletal_muscle" = c("skeletal muscle"), 
              "bone" = c("bone", "femur", "patella", "tibia", "fibula", 
                         "clavicle", "scapula", "humeral", "flexor", 
                         "supinator", "radius", "ulna", "humerus", 
                         "carpus", "metacarpus", "phalanges", "marrow"), 
              "marrow" = c("marrow"), "knee" = c("knee"), "head" = c("head"), 
              "leg" = c("leg"), "arm" = c("arm"), 
              "thorax" = c("thorax"), "chest" = c("chest"),
              "spine" = c("spine", "spinal"), "foot" = c("foot"), 
              "hand" = c("hand"),
              "skin" = c("cutaneous", "skin", "melanocyte"), 
              "kidney" = c("kidney", "renal", "abdominal gland", "nephron"), 
              "corpuscule" = c("corpuscule"), "tubule" = c("tubule"), 
              "gallbladder" = c("gallbladder", "gall bladder", "gall-bladder"),
              "saliva" = c("saliva", "sputum"), "sputum" = c("sputum"), 
              "mucus" = c("mucus"), "mucinous" = c("mucinous"), 
              "fiber" = c("fiber", "fibrous"), "cartilage" = c("cartilage"),
              "joint" = c("joint"),
              "heart" = c("cardiac", "heart", "superior vena cava", "aorta", 
                          "pulmonary artery", "pulmonary vein", "atrium", 
                          "pulmonary valve", "ticuspid valve", 
                          "inferior vena cava", "mitral valve", 
                          "aortic valve", "ventricle"), 
              "esophagus" = c("esophag", "oesophag"), 
              "stomach" = c("stomach", "gastric"), 
              "barretts" = c("( |^)BE( |$|:)", "barretts", "barrett's"), 
              "dysplasia" = c("( |^)(H|L)GD( |$|:)", "dysplasia"),
              "low_grade" = c("low grade"), "high_grade" = c("high_grade"),
              "squamous" = c("( |^)SQ( |$|:)", "squamous"), 
              "colorectal" = c("colorec"), "intestine" = c("colorec"),
              "colon" = c("colon", "colorec", "large intestine", "cecum"), 
              "intestine" = c("colon", "colorec", "large intestine", "cecum"),
              "rectum" = c("colorec", "rectal", "rectum", "anus"),
              "respiratory_system" = c("lung", "bronchi", "alveol", 
                                       "interstiti", "pleura", 
                                       "trachea", "windpipe", "wind pipe", 
                                       "bronchi", "airway"), 
              "lung" = c("lung", "bronchi", "alveol", "interstiti", "pleura"), 
              "alveolar" = c("alveolar"), "lepidic" = c("lepidic"),
              "windpipe" = c("trachea", "windpipe", "wind pipe", "bronchi", 
                             "airway"), 
              "nervous_system" = c("astrocyte", "oligodendrocyte", "ependymal", 
                                   "schwann", "satellite cell", "glia"), 
              "liver" = c("liver", "hepato", "kupff"), 
              "bladder" = c("bladder", "urothel"), 
              "brain" = c("brain", "cerebrum", "cerebral", "cerebellum", 
                          "cerebelli", "dorsolat", "medulla", "lobe", 
                          "prefront", "occipital", "falx", "meningeal", 
                          "supratentorial", "fossa", "sellar", "grey matter", 
                          "gray matter", "white matter", "tentorium", 
                          "tentorial", "cortex", "hippocampus"), 
              "frontal_lobe" = c("frontal lobe"), 
              "frontal_cortex" = c("frontal cortex"), 
              "parietal_lobe" =c("parietal lobe"), 
              "occipital_lobe" = c("occipital lobe"), 
              "prefrontal_lobe" = c("prefront"), 
              "brainstem" = c("brain stem", "brainstem"), 
              "hippocampus" =c("hippocampus"),
              "cerebellum" = c("cerebellum"), 
              "cerebrum" = c("cerebrum", "cerebral"), 
              "cortex" = c("cortex"), "white_matter" = c("white matter"), 
              "gray_matter" = c("gray matter", "grey matter"), 
              "cortex" = c("cortex"), 
              "occipital" = c("occipital"), 
              "placenta" = c("chorion", "villus", "placent"), 
              "umbilical_cord" = c("umbilical", "cord blood"), 
              "uterus" = c("uterus", "uteri", "endometr"), 
              "ovary" = c("ovary", "ovari", "endometrium", "endometrioid"), 
              "fallopian_tube" = c("fallop"), "prostate" = c("prostate"), 
              "neck" = c("neck", "thyroid"), "thyroid_gland" = c("thyroid"), 
              "adrenal" = c("adrenal"), 
              "eye" = c("eye", "uvea", "optic nerve", "cone", "rod", "retina"), 
              "endocrine_system" = c("endocrine", "pineal", "pituitary", 
                                     "pancreas", "pancreat", "adren", 
                                     "thyroid", "hypothalamus", 
                                     "adrenal cortex", "adreno", 
                                     "paraganglioma", "paraganglioma", 
                                     "pheochromocytoma", "zona", 
                                     "glomerulosa", "fasciculata", 
                                     "reticularis", "ovary", "ovari",
                                     "testic", "teste"),
              "pancreas" = c("pancreas", "pancreat"),
              "skin" = c("skin", "epidermis", "keratinocyt"), 
              "keratinocyte" = c("keratinocyt"), "breast" = c("breast"),
              "lymphatic_system" = c("lymph", "spleen", "thymus"), 
              "oral" = c("mouth", "buccal", "labial", "lip", "tongue", 
                         "lingual", "throat", "masticatory"), 
              "throat" = c("throat"), 
              "buccal" = c("buccal", "cheek swab", "mouth swab"), 
              "neuron" = c("neuro", "neural", "nerve", "dendrite", "axon"), 
              "glia" = c("glia"), "epithelial" = c("epithel"),
              "endothelium" = c("endothel"), 
              "stem_cell" = c("stem cell", "pluripot", "ipsc", "iPSC"), 
              "induced_pluripotent_stem_cell" = c("ipsc", "iPSC"), 
              "fibroblast" = c("fibroblast"), 
              "crypt" = c("crypt"),"ectoderm" = c("ectoderm"), 
              "mucosa" = c("mucosa"), "naive" = c("naive", "naïve"), 
              "primed" = c("primed"), 
              "nasal" = c("nasal", "nose", "septum", "sinus"),
              "sperm" = c("sperm", "semen"), "gamete" = c("gamete"),
              "adipose" = c("adipose", "fat", "visceral")); return(txl)
}

# storage condition terms
#' Storage condition terms for storageinfo variable
#' 
#' List containing search terms and string queries for regex pattern matching.
#' 
#' @return List of terms (names) and string queries (values).
#' @export
md_post_storage <- function(){
  sll <- list("frozen" = c("FF$", "frozen", "frzn", "fzn"), 
              "fresh_frozen" = c("FF$", "frozen", "frzn", "fzn"), 
              "FF" = c("FF$", "frozen", "frzn", "fzn"), 
              "formalin_fixed_paraffin_embedded" = c("FFPE", "formalin"), 
              "FFPE" = c("FFPE", "formalin")); return(sll)
}

#' Age terms to seed regex queries
#' 
#' Age terms to seed regex queries, called by md_post_handle_age().
#' 
#' @return List of terms (names) and string queries (values).
#' @export
md_age_infol <- function(){
  ageinfol <- list("adult" = c("adult", "old", "senior"), 
                   "young" = c("neonatal", "pediatric", "prepubescent", 
                               "youth","young","child","infant","postnate"),
                   "prenatal"=c("embryo","embryonic","prenatal","prenate"),
                   "fetal" = c("fetal", "foetal", "fetus"),
                   "neonatal" = c("neonate", "neonatal"),
                   "maternal" = c("maternal", "mother"));return(ageinfol)
}

#' Age unit terms to seed regex queries
#' 
#' Age unit terms to seed regex queries, called by md_post_handle_age().
#' 
#' @return List of unit term labels (names) and string query seed terms 
#' (values).
#' @export
md_age_unitl <- function(){
  ageunitl <- list("years" = c("year", "yr", "y ", "(0-9)y.*"), 
                   "months" = c("month", "mo"), "weeks" = c("week", "wk"), 
                   "days" = c("day", "dy"), "passage" = c("passage"))
  return(ageunitl)
}

# age terms and logic handling
#' Handle age term mappings for metadata postprocessing
#' 
#' This function is called by md_postprocess() in order to handle the logic of
#' age term mappings from preprocessed metadata.
#' 
#' @param mdpre Table of preprocessed metadata.
#' @param mdpost Table of postprocessed metadata.
#' @param mdpost.vl List mapping term categories (names) to column
#' names in the mdpost postprocessed metadata matrix to be generated.
#' @param mdpre.vl List mapping term categories (names) to column
#' names in the mdpre preprocessed metadata matrix.
#' @param verbose Whether to show status messages (TRUE).
#' @seealso md_postprocess(); md_preprocess()
#' @return Postprocessed metadata table with a new column of mapped age
#' info, values.
#' @export
md_post_handle_age <- function(mdpre, mdpost, mdpost.vl = list("age" = "age"),
                               mdpre.vl = list("sample_id" = "gsm",
                                               "sample_title" = "gsm_title",
                                               "info" = "info","age" = "age",
                                               "age_temp" = "age_temp"),
                               verbose = TRUE){
  ap <- mdpre[,mdpre.vl[["age"]]]
  if(verbose){message("Getting filtered, formatted ages...")}
  av <- unlist(lapply(ap, function(x){
    xsplit.num.form <- "NA"
    if(!x == "NA"){
      xval <- unlist(strsplit(x, ";")) # split values
      xvf <- xval[grepl("[0-9]", xval)][1] # catch first numeric value
      xsplit = unlist(strsplit(xvf, "[a-zA-Z]+"))
      xsplit.num <- xsplit[grepl("[0-9]+", xsplit)][1] # catch 1st numeric
      symv <- "[a-zA-Z]| |:|\\)|\\(|/|!|?|\\_|+|," # replace remaining symbols
      xsplit.num.form <- gsub(symv, "", xsplit.num)};return(xsplit.num.form)
  })); av[av == ""] <- "NA"
  mdpost[,mdpost.vl[["age"]]] <- paste0("age_val:", as.character(av))
  if(verbose){message("Getting filtered age data as new variable...")}
  age.val.filt <- unlist(lapply(ap, function(x){
    xvf <- "NA"
    if(!x == "NA"){
      xval <- unlist(strsplit(x, ";")); cond1 <- grepl("[0-9]", xval) 
      cond2 <- grepl("age", xval) & !grepl("stage", xval) # catch age tag
      xvf <- paste(xval[(cond1|cond2)], collapse = ";")};return(xvf)}))
  agetemp.cname <- mdpre.vl[["age_temp"]]; mdpre[,ncol(mdpre) + 1] <- "NA"
  colnames(mdpre)[ncol(mdpre)] <- agetemp.cname
  if(!is.null(age.val.filt)){mdpre[,agetemp.cname] <- age.val.filt}
  if(verbose){message("Adding age metadata...")};ageunitl <- md_age_unitl()
  gf.run <- rep(FALSE, nrow(mdpre))
  for(term in names(ageunitl)){ # allows first units match
    gf.rep <- get_filt(v = get_pstr(v = ageunitl[[term]]), m = mdpre,
                       varl = agetemp.cname)
    mdpost[,mdpost.vl[["age"]]] <- appendvar(var=mdpost.vl[["age"]],
                                                 val=paste0("age_units:",term),
                                                  filtv=gf.rep & !gf.run,
                                                  m = mdpost)
    gf.run <- gf.run|gf.rep};if(verbose){message("Adding age group info...")}
  ageinfol <- md_age_infol(); lgf = list()
  age.cname <- mdpre.vl[["age"]]; title.cname <- mdpre.vl[["sample_title"]]
  which.var <- c(mdpre.vl[["age"]], mdpre.vl[["sample_title"]])
  for(term in names(ageinfol)){
    age.var <- get_pstr(v = ageinfol[[term]])
    lgf[[term]][[age.cname]]<-get_filt(v=age.var,m=mdpre,varl=which.var)
    title.var <- get_pstr(v = ageinfol[[term]])
    lgf[[term]][[title.cname]]<-get_filt(v=title.var,m=mdpre,varl=which.var)}
  if(verbose){message("Handling age info map logic...")};termv<-names(ageinfol)
  for(r in seq(nrow(mdpost))){
    sv.term <- "NA"; bool.term.age <- bool.term.title <- c()
    for(t in termv){
      bool.term.age <- c(bool.term.age, lgf[[t]][[age.cname]][r])
      bool.term.title <- c(bool.term.title, lgf[[t]][[title.cname]][r])}
    tf.age<-termv[which(bool.term.age)];tf.title<-termv[which(bool.term.title)]
    sv.term <- ifelse(length(tf.age) == 1, tf.age, 
                      ifelse(length(tf.title) == 1, tf.title, "NA"))
    info.val <- paste0("age_info:", sv.term)
    mdpost[,mdpost.vl[["age"]]][r] <- paste0(mdpost[,mdpost.vl[["age"]]][r],
                                             ";", info.val)}
  return(mdpost)
}

# main postprocess function

#' Postprocess metadata prepared using md_preprocess()
#'
#' Perform postprocessing of previously prepreocessed sample metadata. This 
#' produces the new harmonized variables (specified by arg mdpost.vl),
#' where harmonization means terms are mapped, lowercase, and "_" separated.
#' Variable entries can include multiple terms separated by ";". The args
#' mdpre.vl and mdpost.vl specify the various variable titles in
#' the preprocess and postprocess dataset. The vars disease.search.vars,
#' tissue.search.vars, and storage.info.vars specify the mdpre variables to
#' search for disease, tissue, and storage info mappings.
#' 
#' @param ts The timestamp for this run.
#' @param mdpre The matrix containing preprocessed metadata (returned from
#' md_preprocess()).
#' @param mdpost.fname Filename for newly mapped postprocessed metadata.
#' @param md.dpath Path to the directory containing the preprocessed 
#' metadata matrix, where newly postprocessed metadata will be stored.
#' @param mdpre.vl List mapping term categories (names) to column
#' names in the mdpre preprocessed metadata matrix.
#' @param mdpost.vl List mapping term categories (names) to column
#' names in the mdpost postprocessed metadata matrix to be generated.
#' @param disease.search.vars Term categories to search in preprocessed 
#' metadata for the disease term mappings
#' @param tissue.search.vars Term categories in preprocessed metadata to
#' search for tissue term mappings.
#' @param storage.info.vars Term categories in preprocessed metadata to
#' search for storage information term mappings.
#' @param verbose Whether to show status messages (TRUE).
#' @seealso md_preprocess()
#' @return Postprocessed metadata table.
#' @export
md_postprocess <- function(ts, mdpre, mdpost.fname = "md_postprocess",
                           md.dpath = file.path("recount-methylation-files", 
                                              "metadata"),
                           mdpre.vl=list("study_id"="gse","sample_id"="gsm",
                                         "sample_title" = "gsm_title",
                                         "disease" = "disease_state",
                                         "sample_type" = "sample_type",
                                         "sex" = "sex", "info" = "info",
                                         "age" = "age", 
                                         "age_temp" = "age_temp"),
                           mdpost.vl = list("tissue" = "tissue",
                                            "disease" = "disease",
                                            "age" = "age", "sex" = "sex",
                                            "storageinfo" = "storageinfo"),
                           disease.search.vars = c("sample_title", "disease"),
                           tissue.search.vars = c("sample_type", 
                                                  "sample_title"),
                           storage.info.vars = c("sample_type", 
                                                 "sample_title", "info"),
                           verbose = TRUE){
  mdpost.fn <- paste0(mdpost.fname, "_", ts, ".rda")
  mdpost.fpath <- file.path(md.dpath, mdpost.fn)
  if(verbose){message("Will save mdpost data to ", mdpost.fpath, "...")}
  mdpost <- mdpre[,c(mdpre.vl[["sample_id"]], mdpre.vl[["study_id"]], 
                     mdpre.vl[["sample_title"]])]
  mdpost[,mdpost.vl[["tissue"]]] <- mdpost[,mdpost.vl[["disease"]]] <- "NA"
  mdpost[,mdpost.vl[["age"]]] <- mdpost[,mdpost.vl[["sex"]]] <- "NA"
  mdpost[,mdpost.vl[["storageinfo"]]] <- "NA"
  if(verbose){message("Getting disease status...")}; dxl <- md_post_disease()
  which.var <- unlist(mdpre.vl[names(mdpre.vl) %in% disease.search.vars])
  for(dx in names(dxl)){
    ssv <- dxl[[dx]]; pstr <- get_pstr(v = ssv)
    gfilt <- get_filt(v = pstr, m = mdpre, ntfilt = ssv, varl = which.var)
    dx.var <- appendvar(var = mdpost.vl[["disease"]], val = dx, 
                        filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["disease"]]] <- dx.var}
  if(verbose){message("Getting disease terms for cancers by type/location...")}
  which.var <- unlist(mdpre.vl[names(mdpre.vl) %in% disease.search.vars])
  cxsubl <- md_post_cancertype()
  for(cxsub in names(cxsubl)){
    ssv <- cxsubl[[cxsub]]; pstr <- get_pstr(v = ssv)
    gfilt <- get_filt(v = pstr, m = mdpre, varl = which.var)
    dxvar1 <- appendvar(var = mdpost.vl[["disease"]], val = cxsub, 
                        filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["disease"]]] <- dxvar1
    dxvar2 <- appendvar(var = mdpost.vl[["disease"]], val = "cancer", 
                        filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["disease"]]] <- dxvar2}
  if(verbose){message("Getting leukemia disease terms...")};leukl <- md_post_leukemia()
  which.var <- unlist(mdpre.vl[names(mdpre.vl) %in% disease.search.vars])
  for(leuk in names(leukl)){
    pstr <- suppressMessages(get_pstr(v = leukl[[leuk]]))
    gfilt <- suppressMessages(get_filt(v = pstr, m = mdpre, ntfilt = pstr,
                                     varl = which.var))
    dxvar1 <- appendvar(var = mdpost.vl[["disease"]], val = leuk, 
                        filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["disease"]]] <- dxvar1
    dxvar2 <- appendvar(var = mdpost.vl[["disease"]], val = "cancer", 
                        filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["disease"]]] <- dxvar2}
  if(verbose){message("Getting tissue and disease for cancer subtypes...")}
  which.var <- c(mdpre.vl[["sample_title"]], mdpre.vl[["sample_type"]],
                 mdpre.vl[["disease"]]);cxl <- md_post_cancer()
  for(cx in names(cxl)){
    ssv <- cxl[[cx]]; pstr <- get_pstr(v = ssv)
    gfilt<-suppressMessages(get_filt(v=pstr,m=mdpre,ntfilt=ssv,varl=which.var))
    txvar <- appendvar(var = mdpost.vl[["tissue"]], val = cx, 
                       filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["tissue"]]] <- txvar
    dxvar <- appendvar(var = mdpost.vl[["disease"]], val = "cancer", 
                       filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["disease"]]] <- dxvar}
  if(verbose){message("Getting tissue terms for cancers by type/location...")}
  which.var <- unlist(mdpre.vl[names(mdpre.vl) %in% tissue.search.vars])
  for(cxsub in names(cxsubl)){
    ssv <- cxsubl[[cxsub]]; pstr <- get_pstr(v = ssv)
    gfilt <- get_filt(v = pstr, m = mdpre, varl = which.var);
    txvar1 <- appendvar(var = mdpost.vl[["tissue"]], val = cxsub, 
                        filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["tissue"]]] <- txvar1
    txvar2 <- appendvar(var = mdpost.vl[["tissue"]], val = "cancer", 
                        filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["tissue"]]] <- txvar2}
  if(verbose){message("Getting leukemia tissue terms...")}
  which.var <- unlist(mdpre.vl[names(mdpre.vl) %in% tissue.search.vars])
  for(leuk in names(leukl)){
    pstr <- get_pstr(v = leukl[[leuk]])
    gfilt <- get_filt(v = pstr, m = mdpre, varl = which.var)
    txvar <- appendvar(var = mdpost.vl[["tissue"]], val = leuk, 
                       filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["tissue"]]] <- txvar}
  if(verbose){message("Getting tissue annotations...")}
  which.var <- unlist(mdpre.vl[names(mdpre.vl) %in% tissue.search.vars])
  txl <- list();txl[["blood"]] <- md_post_tissue_blood()
  txl[["other"]] <- md_post_tissue()
  for(sublist in txl){
    for(tx in names(sublist)){
      pstr <- get_pstr(v = sublist[[tx]])
      gfilt <- get_filt(v = pstr, m = mdpre, varl = which.var)
      txvar <- appendvar(var = mdpost.vl[["tissue"]], val = tx, 
                         filtv = gfilt, m = mdpost)
      mdpost[,mdpost.vl[["tissue"]]] <- txvar}}
  if(verbose){message("Getting storage info...")}
  which.var <- unlist(mdpre.vl[names(mdpre.vl) %in% storage.info.vars])
  lstorage <- md_post_storage()
  for(sn in names(lstorage)){
    pstr <- get_pstr(v = lstorage[[sn]])
    gfilt <- get_filt(v = pstr, m = mdpre, varl = which.var)
    sinfovar <- appendvar(var = mdpost.vl[["storageinfo"]], val = sn, 
                          filtv = gfilt, m = mdpost)
    mdpost[,mdpost.vl[["storageinfo"]]] <- sinfovar}
  if(verbose){message("Getting age info...")}
  mdpost<-suppressMessages(md_post_handle_age(mdpre=mdpre,mdpost=mdpost,
                                              mdpre.vl=mdpre.vl,
                                              mdpost.vl = mdpost.vl,
                                              verbose = verbose))
  if(verbose){message("Getting sex info...")}
  mdpre.sex.cname <- mdpre.vl[["sex"]];mdpost.sex.cname <- mdpost.vl[["sex"]]
  femv <- get_pstr(v = c("female", "f", "FEMALE"))
  malev <- get_pstr(v = c("male", "MALE", "m"))
  mdpost[,mdpost.sex.cname] <- ifelse(grepl(femv, mdpre[,mdpre.sex.cname]),"F",
                       ifelse(grepl(malev, mdpre[,mdpre.sex.cname]),"M", "NA"))
  if(verbose){message("Saving mdpost to ", mdpost.fpath)}
  save(mdpost, file = mdpost.fpath);return(NULL)
}

