library(cBioPortalData)

cbio <- cBioPortal()
studies <- getStudies(cbio, buildReport = TRUE)

setCache(here::here("data/cbioportal_cache/"))

tcga_hnsc_tmb <- clinicalData(cbio, studyId = "hnsc_tcga") |> 
    dplyr::select(patient = patientId,
                  sample = sampleId,
                  tmb = TMB_NONSYNONYMOUS) |> 
    dplyr::mutate(tmb = as.numeric(tmb))
tcga_hnsc_pancan_tmb <- clinicalData(
    cbio, 
    studyId = "hnsc_tcga_pan_can_atlas_2018") |> 
    dplyr::select(patient = patientId,
                  sample = sampleId,
                  tmb_pancan = TMB_NONSYNONYMOUS) |> 
    dplyr::mutate(tmb_pancan = as.numeric(tmb_pancan))
tcga_hnsc_tmb <- dplyr::full_join(
    tcga_hnsc_tmb, tcga_hnsc_pancan_tmb,
    join_by(sample == sample, patient == patient)
)

saveRDS(tcga_hnsc_tmb, here::here("output/data/tcga_hnsc_tmb.rds"))
