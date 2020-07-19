library(rmongodb)
library(mongolite)
library(future)
library(dplyr)

# Cargamos la Coleccion HGDP de MongoDB
HGDP_collection <- mongo(collection = "human", db = "HGDP")


# CARGANDO FUNCIONES
search_mongo <- function(IDs){
  a <- c('"V1"', IDs)
  a <- paste0(a, " : 1", collapse = ",")
  search <- paste0('{', a, '}')
  df <- HGDP_collection$find('{}', fields = search)
  return(df)

}

primer_limpia <- function(x){
  c <- list()
  for(i in 1:ncol(x)){
    vect <- paste(x[,i], collapse = "")
    vect <- stringr::str_remove_all(vect, "-")
    names(vect) <- colnames(x)[i]
    c[[i]] <- vect

  }

  return(c)
}

segunda_limpia <- function(letters_n, total_pop){

  while(letters_n != ""){
    t <- stringr::str_count(letters_n, "T")
    c <- stringr::str_count(letters_n, "C")
    g <- stringr::str_count(letters_n, "G")
    a <- stringr::str_count(letters_n, "A")

    all_nucl <- c(t,c,g,a)
    names(all_nucl) <- c("T", "C", "G", "A")
    all_nucl <- all_nucl[all_nucl > 0]

    if(length(all_nucl) == 1 ){
      freq_ancestral <- (all_nucl[which.max(all_nucl)] / 2) / total_pop
      freq_minor <- 0
      ancestral <- names(freq_ancestral)
      minor <- names(freq_ancestral)
      freq_ancestral <- unname(freq_ancestral)
      # freq_minor <- unname(freq_minor)
      data_freqs <- data.frame(SNP = names(letters_n),
                               Ancestral_freq = freq_ancestral,
                               MAF = 0,
                               Ancestral_allel = ancestral,
                               Minor_allel = minor)

    } else {

      freq_ancestral <- (all_nucl[which.max(all_nucl)] / 2) / total_pop
      freq_minor <- (all_nucl[which.min(all_nucl)] / 2) / total_pop
      ancestral <- names(freq_ancestral)
      minor <- names(freq_minor)
      freq_ancestral <- unname(freq_ancestral)
      freq_minor <- unname(freq_minor)
      data_freqs <- data.frame(SNP = names(letters_n),
                               Ancestral_freq = freq_ancestral,
                               MAF = freq_minor,
                               Ancestral_allel = ancestral,
                               Minor_allel = minor)
    }
    return(data_freqs)

  }
}

tercer_limpia <- function(contents_request){

  contents_request[sapply(contents_request, is.null)] <- NULL
  prueba4 <- rlist::list.stack(contents_request)
  return(prueba4)

}

limpieza_total <- function(df){
  total_pop <- ncol(df) -2
  prueba1 <- t(df[,-1])
  colnames(prueba1) <- (prueba1[1,])
  prueba1 <- prueba1[-1,]
  prueba2 <- primer_limpia(prueba1)
  names(prueba2) <- colnames(prueba1)

  future::plan(multiprocess)
  content <- furrr::future_map(prueba2, purrr::safely(segunda_limpia),
                               total_pop = total_pop,
                               .progress = TRUE)
  contents_1 <- purrr::transpose(content)
  contents_request <- contents_1[["result"]]
  prueba3 <- tercer_limpia(contents_request)

  return(prueba3)
}


# Adygei in Caucasus, Russia HGDP

IDs_adygei <- c('"HGDP01403"', '"HGDP01383"', '"HGDP01388"', '"HGDP01404"', '"HGDP01384"','"HGDP01396"',
                '"HGDP01385"', '"HGDP01397"', '"HGDP01400"', '"HGDP01381"', '"HGDP01386"', '"HGDP01398"',
                '"HGDP01382"', '"HGDP01387"', '"HGDP01399"')


Adygei_russia <- search_mongo(IDs_adygei) %>%
  limpieza_total()

usethis::use_data(Adygei_russia, overwrite = TRUE)

# Balochi in Pakistan HGDP

IDs_balochi <- c('"HGDP00056"', '"HGDP00068"', '"HGDP00070"', '"HGDP00082"', '"HGDP00094"', '"HGDP00052"',
                 '"HGDP00057"', '"HGDP00076"', '"HGDP00088"', '"HGDP00064"', '"HGDP00060"', '"HGDP00072"',
                 '"HGDP00096"', '"HGDP00078"', '"HGDP00080"', '"HGDP00092"', '"HGDP00054"', '"HGDP00066"',
                 '"HGDP00086"', '"HGDP00098"', '"HGDP00062"', '"HGDP00074"')
Balochi_Pakistan <- search_mongo(IDs_balochi) %>%
  limpieza_total()

usethis::use_data(Balochi_Pakistan, overwrite = TRUE)

# -------
# Bantu in kenya HGDP

IDs_bandu <- c('"HGDP01408"', '"HGDP01415"', '"HGDP01411"', '"HGDP01416"', '"HGDP01405"', '"HGDP01412"',
               '"HGDP01418"', '"HGDP01406"', '"HGDP01413"', '"HGDP01419"')
Bandu_Kenya <- search_mongo(IDs_bandu) %>%
  limpieza_total()

usethis::use_data(Bandu_Kenya, overwrite = TRUE)


# Bandu in South Africa HGDP

IDs_bandu_southafrica <- c('"HGDP00993"', '"HGDP01031"', '"HGDP00994"', '"HGDP01033"')
Bandu_SouthAfrica <- search_mongo(IDs_bandu_southafrica) %>%
  limpieza_total()

usethis::use_data(Bandu_SouthAfrica, overwrite = TRUE)


# Basque in France HGDP

IDs_basque <- c('"HGDP01357"', '"HGDP01369"', '"HGDP01371"', '"HGDP01376"', '"HGDP01358"', '"HGDP01360"',
                '"HGDP01372"', '"HGDP01377"', '"HGDP01359"', '"HGDP01361"', '"HGDP01378"', '"HGDP01380"',
                '"HGDP01366"', '"HGDP01373"', '"HGDP01374"', '"HGDP01379"', '"HGDP01362"', '"HGDP01367"',
                '"HGDP01375"', '"HGDP01363"', '"HGDP01368"', '"HGDP01370"')
Basque_France <- search_mongo(IDs_basque) %>%
  limpieza_total()

usethis::use_data(Basque_France, overwrite = TRUE)



# Bedouin in Negev, Israel HGDP

IDs_bedouin <- c('"HGDP00607"', '"HGDP00614"', '"HGDP00619"', '"HGDP00640"', '"HGDP00645"', '"HGDP00621"',
                 '"HGDP00626"', '"HGDP00638"', '"HGDP00608"', '"HGDP00610"', '"HGDP00615"', '"HGDP00634"',
                 '"HGDP00639"', '"HGDP00641"', '"HGDP00646"', '"HGDP00622"', '"HGDP00627"', '"HGDP00653"',
                 '"HGDP00630"', '"HGDP00635"', '"HGDP00623"', '"HGDP00628"', '"HGDP00609"', '"HGDP00611"',
                 '"HGDP00654"', '"HGDP00642"', '"HGDP00647"', '"HGDP00643"', '"HGDP00648"', '"HGDP00631"',
                 '"HGDP00636"', '"HGDP00624"', '"HGDP00629"', '"HGDP00612"', '"HGDP00701"', '"HGDP00649"',
                 '"HGDP00651"', '"HGDP00613"', '"HGDP00618"', '"HGDP00620"', '"HGDP00632"', '"HGDP00637"',
                 '"HGDP00644"', '"HGDP00625"')
Bedouin_Israel <- search_mongo(IDs_bedouin) %>%
  limpieza_total()

usethis::use_data(Bedouin_Israel, overwrite = TRUE)


# Bergamo Italian in Bergamo, Italy HGDP

IDs_bergamo <- c('"HGDP01174"', '"HGDP01155"', '"HGDP01149"', '"HGDP01151"', '"HGDP01156"', '"HGDP01152"',
                 '"HGDP01157"', '"HGDP01171"', '"HGDP01177"', '"HGDP01173"')
Bergamo_Italy <- search_mongo(IDs_bergamo) %>%
  limpieza_total()

usethis::use_data(Bergamo_Italy, overwrite = TRUE)


# Biaka in Central African Republic HGDP

IDs_biaka <- c('"HGDP00455"', '"HGDP00479"', '"HGDP00981"', '"HGDP00986"', '"HGDP01092"', '"HGDP00470"',
               '"HGDP00475"', '"HGDP01086"', '"HGDP00464"', '"HGDP00469"', '"HGDP00452"', '"HGDP01094"',
               '"HGDP01087"', '"HGDP00453"', '"HGDP00458"', '"HGDP00460"', '"HGDP00465"', '"HGDP00472"',
               '"HGDP01090"', '"HGDP00454"', '"HGDP00459"', '"HGDP00466"', '"HGDP00473"', '"HGDP00985"')
Bianka_CentralAfricaRepublic <- search_mongo(IDs_biaka) %>%
  limpieza_total()

usethis::use_data(Bianka_CentralAfricaRepublic, overwrite = TRUE)


# Bougainville in Bougainville HGDP

IDs_bougainville <- c('"HGDP00664"', '"HGDP00823"', '"HGDP00490"', '"HGDP00661"', '"HGDP00491"', '"HGDP00655"',
                      '"HGDP00662"', '"HGDP00787"', '"HGDP00663"', '"HGDP00788"', '"HGDP01027"')
Bougainville_Bougainville <- search_mongo(IDs_bougainville) %>%
  limpieza_total()

usethis::use_data(Bougainville_Bougainville, overwrite = TRUE)


# Brahui in Pakistan HGDP

IDs_brahui <- c('"HGDP00001"', '"HGDP00049"', '"HGDP00013"', '"HGDP00025"', '"HGDP00037"', '"HGDP00007"',
                '"HGDP00045"', '"HGDP00033"', '"HGDP00021"', '"HGDP00003"', '"HGDP00039"', '"HGDP00041"',
                '"HGDP00015"', '"HGDP00009"', '"HGDP00011"', '"HGDP00023"', '"HGDP00047"', '"HGDP00035"',
                '"HGDP00029"', '"HGDP00031"', '"HGDP00043"', '"HGDP00005"', '"HGDP00017"')
Brahui_Pakistan <- search_mongo(IDs_brahui) %>%
  limpieza_total()

usethis::use_data(Brahui_Pakistan, overwrite = TRUE)


# Burusho in Pakistan HGDP

IDs_burusho <- c('"HGDP00392"', '"HGDP00397"', '"HGDP00412"', '"HGDP00417"', '"HGDP00359"', '"HGDP00444"',
                 '"HGDP00402"', '"HGDP00407"', '"HGDP00433"', '"HGDP00438"', '"HGDP00445"', '"HGDP00351"',
                 '"HGDP00356"', '"HGDP00382"', '"HGDP00364"', '"HGDP00371"', '"HGDP00376"', '"HGDP00388"',
                 '"HGDP00341"', '"HGDP00346"', '"HGDP00372"', '"HGDP00423"')
Burucho_Pakistan <- search_mongo(IDs_burusho) %>%
  limpieza_total()

usethis::use_data(Burucho_Pakistan, overwrite = TRUE)

# -----------
# Cambodian in Cambodia HGDP

IDs_cambodia <- c('"HGDP00715"', '"HGDP00711"', '"HGDP00716"', '"HGDP00712"', '"HGDP00720"', '"HGDP00714"',
                  '"HGDP00719"', '"HGDP00721"')
Cambodian_Cambodia <- search_mongo(IDs_cambodia) %>%
  limpieza_total()

usethis::use_data(Cambodian_Cambodia, overwrite = TRUE)


# Colombian in Colombia HGDP

IDs_colombian <- c('"HGDP00703"', '"HGDP00708"', '"HGDP00710"', '"HGDP00704"', '"HGDP00970"')
Colombian_Colombia <- search_mongo(IDs_colombian) %>%
  limpieza_total()

usethis::use_data(Colombian_Colombia, overwrite = TRUE)


# Dai in China HGDP

IDs_Dai <- c('"HGDP01310"', '"HGDP01309"', '"HGDP01311"', '"HGDP01316"', '"HGDP01313"')
Dai_China <- search_mongo(IDs_Dai) %>%
  limpieza_total()

usethis::use_data(Dai_China, overwrite = TRUE)


# Daur in China HGDP

IDs_daur <- c('"HGDP01213"', '"HGDP01218"', '"HGDP01220"', '"HGDP01214"', '"HGDP01219"', '"HGDP01221"',
              '"HGDP01222"', '"HGDP01216"', '"HGDP01217"')
Daur_China <- search_mongo(IDs_daur) %>%
  limpieza_total()

usethis::use_data(Daur_China, overwrite = TRUE)


# Druze in Carmel, Israel HGDP

IDs_druze <- c('"HGDP00563"', '"HGDP00599"', '"HGDP00602"', '"HGDP00582"', '"HGDP00587"', '"HGDP00594"',
               '"HGDP00575"', '"HGDP00568"', '"HGDP00583"', '"HGDP00588"', '"HGDP00590"', '"HGDP00595"',
               '"HGDP00564"', '"HGDP00571"', '"HGDP00557"', '"HGDP00576"', '"HGDP00558"', '"HGDP00560"',
               '"HGDP00565"', '"HGDP00572"', '"HGDP00584"', '"HGDP00591"', '"HGDP00577"', '"HGDP00604"',
               '"HGDP00578"', '"HGDP00580"', '"HGDP00559"', '"HGDP00561"', '"HGDP00566"', '"HGDP00573"',
               '"HGDP00600"', '"HGDP00562"', '"HGDP00579"', '"HGDP00581"', '"HGDP00567"', '"HGDP00574"',
               '"HGDP00598"', '"HGDP00601"', '"HGDP00606"', '"HGDP00586"')
Druze_Israel <- search_mongo(IDs_druze) %>%
  limpieza_total()

usethis::use_data(Druze_Israel, overwrite = TRUE)

# -----------
# French in France HGDP

IDs_french <- c('"HGDP00513"', '"HGDP00518"', '"HGDP00520"', '"HGDP00525"', '"HGDP00537"', '"HGDP00538"',
                '"HGDP00514"', '"HGDP00519"', '"HGDP00515"', '"HGDP00534"', '"HGDP00539"', '"HGDP00522"',
                '"HGDP00527"', '"HGDP00511"', '"HGDP00528"', '"HGDP00535"', '"HGDP00516"', '"HGDP00523"',
                '"HGDP00512"', '"HGDP00529"', '"HGDP00531"', '"HGDP00536"', '"HGDP00517"', '"HGDP00524"')
French_France <- search_mongo(IDs_french) %>%
  limpieza_total()

usethis::use_data(French_France, overwrite = TRUE)


# Han in China HGDP

IDs_han <- c('"HGDP00784"', '"HGDP00777"', '"HGDP00811"', '"HGDP00974"', '"HGDP01023"', '"HGDP00812"',
             '"HGDP00817"', '"HGDP00780"', '"HGDP00975"', '"HGDP01024"', '"HGDP00774"', '"HGDP00779"',
             '"HGDP00813"', '"HGDP00818"', '"HGDP00820"', '"HGDP00781"', '"HGDP00786"', '"HGDP00971"',
             '"HGDP00976"', '"HGDP00782"', '"HGDP00821"', '"HGDP00814"', '"HGDP00819"', '"HGDP01021"',
             '"HGDP00977"', '"HGDP00972"', '"HGDP00776"', '"HGDP00815"', '"HGDP00822"', '"HGDP00973"')
Han_China <- search_mongo(IDs_han) %>%
  limpieza_total()
Han_China <- readRDS("Han_China.rds")
usethis::use_data(Han_China, overwrite = TRUE)


# Hazara in Pakistan HGDP

IDs_hazara <- c('"HGDP00099"', '"HGDP00102"', '"HGDP00119"', '"HGDP00121"', '"HGDP00115"', '"HGDP00122"',
                '"HGDP00127"', '"HGDP00103"', '"HGDP00108"', '"HGDP00110"', '"HGDP00104"', '"HGDP00109"',
                '"HGDP00100"', '"HGDP00105"', '"HGDP00129"', '"HGDP00118"', '"HGDP00120"', '"HGDP00106"')
Hazara_Pakistan <- search_mongo(IDs_hazara) %>%
  limpieza_total()

usethis::use_data(Hazara_Pakistan, overwrite = TRUE)


# Hezhen in China HGDP

IDs_hezhen <- c('"HGDP01237"', '"HGDP01233"', '"HGDP01238"', '"HGDP01234"', '"HGDP01239"', '"HGDP01241"',
                '"HGDP01236"')
Hezhen_China <- search_mongo(IDs_hezhen) %>%
  limpieza_total()

usethis::use_data(Hezhen_China, overwrite = TRUE)


# Japanese in Japan HGDP

IDs_japanese <- c('"HGDP00753"', '"HGDP00758"', '"HGDP00760"', '"HGDP00765"', '"HGDP00828"', '"HGDP00791"',
                  '"HGDP00772"', '"HGDP00759"', '"HGDP00761"', '"HGDP00747"', '"HGDP00754"', '"HGDP00766"',
                  '"HGDP00767"', '"HGDP00748"', '"HGDP00750"', '"HGDP00755"', '"HGDP00762"', '"HGDP00751"',
                  '"HGDP00756"', '"HGDP00763"', '"HGDP00768"', '"HGDP00752"', '"HGDP00769"', '"HGDP00771"',
                  '"HGDP00790"', '"HGDP00757"', '"HGDP00764"')
Japanese_Japan <- search_mongo(IDs_japanese) %>%
  limpieza_total()

usethis::use_data(Japanese_Japan, overwrite = TRUE)

# Kalash in Pakistan HGDP

IDs_kalash <- c('"HGDP00304"', '"HGDP00309"', '"HGDP00311"', '"HGDP00277"', '"HGDP00323"', '"HGDP00330"',
                '"HGDP00285"', '"HGDP00274"', '"HGDP00279"', '"HGDP00281"', '"HGDP00298"', '"HGDP00313"',
                '"HGDP00319"', '"HGDP00321"', '"HGDP00326"', '"HGDP00333"', '"HGDP00302"', '"HGDP00307"',
                '"HGDP00315"', '"HGDP00288"', '"HGDP00290"')
Kalash_Pakistan <- search_mongo(IDs_kalash) %>%
  limpieza_total()

usethis::use_data(Kalash_Pakistan, overwrite = TRUE)

# Karitiana in Brazil HGDP

IDs_karitiana <- c('"HGDP01009"', '"HGDP00999"', '"HGDP00995"', '"HGDP01001"', '"HGDP01006"',
                   '"HGDP01013"', '"HGDP01014"', '"HGDP01019"', '"HGDP01010"')
Karitiana_Brazil <- search_mongo(IDs_karitiana) %>%
  limpieza_total()

usethis::use_data(Karitiana_Brazil, overwrite = TRUE)


# Lahu in China HGDP

IDs_lahu <- c('"HGDP01319"', '"HGDP01321"', '"HGDP01326"', '"HGDP01322"', '"HGDP01317"', '"HGDP01318"')
Lahu_China <- search_mongo(IDs_lahu) %>%
  limpieza_total()

usethis::use_data(Lahu_China, overwrite = TRUE)


# Makrani in Pakistan HGDP

IDs_makrani <- c('"HGDP00133"', '"HGDP00140"', '"HGDP00145"', '"HGDP00153"', '"HGDP00158"', '"HGDP00134"',
                 '"HGDP00139"', '"HGDP00141"', '"HGDP00146"', '"HGDP00130"', '"HGDP00135"', '"HGDP00154"',
                 '"HGDP00161"', '"HGDP00136"', '"HGDP00143"', '"HGDP00148"', '"HGDP00150"', '"HGDP00155"',
                 '"HGDP00131"', '"HGDP00137"', '"HGDP00144"', '"HGDP00149"', '"HGDP00151"')
Makrani_Pakistan <- search_mongo(IDs_makrani) %>%
  limpieza_total()

usethis::use_data(Makrani_Pakistan, overwrite = TRUE)



# Mandenka in Senegal HGDP

IDs_mandenka <- c('"HGDP00905"', '"HGDP00912"', '"HGDP00917"', '"HGDP01201"', '"HGDP00906"', '"HGDP00913"',
                  '"HGDP00918"', '"HGDP01202"', '"HGDP01283"', '"HGDP00907"', '"HGDP00914"', '"HGDP00919"',
                  '"HGDP00908"', '"HGDP00910"', '"HGDP01285"', '"HGDP00904"', '"HGDP00909"', '"HGDP00911"',
                  '"HGDP00916"', '"HGDP01200"')
Mandenka_Senegal <- search_mongo(IDs_mandenka) %>%
  limpieza_total()

usethis::use_data(Mandenka_Senegal, overwrite = TRUE)


# Mayan HGDP Mexico

IDs_Maya_HGDP <- c('"HGDP00854"', '"HGDP00859"', '"HGDP00861"', '"HGDP00873"', '"HGDP00862"',
                   '"HGDP00868"', '"HGDP00870"', '"HGDP00875"', '"HGDP00856"', '"HGDP00863"',
                   '"HGDP00869"', '"HGDP00871"', '"HGDP00876"', '"HGDP00864"', '"HGDP00858"',
                   '"HGDP00860"', '"HGDP00865"', '"HGDP00872"', '"HGDP00877"')

Mayan_Mexico <- search_mongo(IDs_Maya_HGDP) %>%
  limpieza_total()

usethis::use_data(Mayan_Mexico, overwrite = TRUE)

# Mbuti in Democratic Republic of Congo HGDP

IDs_mbuti <- c('"HGDP00450"', '"HGDP00462"', '"HGDP00467"', '"HGDP00463"', '"HGDP00468"', '"HGDP01081"',
               '"HGDP00471"', '"HGDP00983"', '"HGDP00984"', '"HGDP00478"')
Mbuti_DemocraticRepublicOfCongo <- search_mongo(IDs_mbuti) %>%
  limpieza_total()

usethis::use_data(Mbuti_DemocraticRepublicOfCongo, overwrite = TRUE)


# Miao in China HGDP

IDs_miao <- c('"HGDP01193"', '"HGDP01194"', '"HGDP01190"', '"HGDP01195"', '"HGDP01189"', '"HGDP01196"',
              '"HGDP01192"', '"HGDP01197"')
Miao_China <- search_mongo(IDs_miao) %>%
  limpieza_total()

usethis::use_data(Miao_China, overwrite = TRUE)


# Mongolian in China HGDP

IDs_mongolian <- c('"HGDP01225"', '"HGDP01232"', '"HGDP01226"', '"HGDP01227"', '"HGDP01230"',
                   '"HGDP01229"', '"HGDP01231"', '"HGDP01224"')
Mongolian_China <- search_mongo(IDs_mongolian) %>%
  limpieza_total()

usethis::use_data(Mongolian_China, overwrite = TRUE)


# Mozabite in Algeria HGDP

IDs_mozabite <- c('"HGDP01256"', '"HGDP01263"', '"HGDP01282"', '"HGDP01268"', '"HGDP01270"', '"HGDP01275"',
                  '"HGDP01264"', '"HGDP01269"', '"HGDP01271"', '"HGDP01276"', '"HGDP01257"', '"HGDP01265"',
                  '"HGDP01272"', '"HGDP01277"', '"HGDP01258"', '"HGDP01260"', '"HGDP01254"', '"HGDP01259"',
                  '"HGDP01261"', '"HGDP01266"', '"HGDP01273"', '"HGDP01280"', '"HGDP01255"', '"HGDP01262"',
                  '"HGDP01267"', '"HGDP01279"')
Mozabite_Algeria <- search_mongo(IDs_mozabite) %>%
  limpieza_total()

usethis::use_data(Mozabite_Algeria, overwrite = TRUE)


# Naxi in China HGDP

IDs_naxi <- c('"HGDP01340"', '"HGDP01341"', '"HGDP01346"', '"HGDP01339"', '"HGDP01342"', '"HGDP01337"')
Naxi_China <- search_mongo(IDs_naxi) %>%
  limpieza_total()

usethis::use_data(Mozabite_Algeria, overwrite = TRUE)


# Northern Han in China HGDP

IDs_northern <- c('"HGDP01287"', '"HGDP01294"', '"HGDP01288"', '"HGDP01290"', '"HGDP01295"', '"HGDP01289"',
                  '"HGDP01291"', '"HGDP01296"', '"HGDP01292"', '"HGDP01293"')
Northern_Han_China <- search_mongo(IDs_northern) %>%
  limpieza_total()

usethis::use_data(Northern_Han_China, overwrite = TRUE)


# Orcadian in Orkney HGDP

IDs_orcadian <- c('"HGDP00804"', '"HGDP00797"', '"HGDP00800"', '"HGDP00805"', '"HGDP00806"', '"HGDP00794"',
                  '"HGDP00799"', '"HGDP00802"', '"HGDP00807"', '"HGDP00795"', '"HGDP00803"', '"HGDP00808"',
                  '"HGDP00810"')
Orcadian_Orkney <- search_mongo(IDs_orcadian) %>%
  limpieza_total()

usethis::use_data(Orcadian_Orkney, overwrite = TRUE)

# Oroqen in China HGDP

IDs_oroqen <- c('"HGDP01206"', '"HGDP01207"', '"HGDP01208"', '"HGDP01204"', '"HGDP01209"', '"HGDP01205"',
                '"HGDP01212"')
Oroqen_China <- search_mongo(IDs_oroqen) %>%
  limpieza_total()

usethis::use_data(Oroqen_China, overwrite = TRUE)


# Palestinian in Israel HGDP

IDs_palestinian <- c('"HGDP00676"', '"HGDP00739"', '"HGDP00741"', '"HGDP00746"', '"HGDP00683"', '"HGDP00688"',
                     '"HGDP00690"', '"HGDP00727"', '"HGDP00734"', '"HGDP00723"', '"HGDP00730"', '"HGDP00735"',
                     '"HGDP00696"', '"HGDP00677"', '"HGDP00684"', '"HGDP00689"', '"HGDP00691"', '"HGDP00678"',
                     '"HGDP00680"', '"HGDP00685"', '"HGDP00692"', '"HGDP00736"', '"HGDP00697"', '"HGDP00700"',
                     '"HGDP00724"', '"HGDP00729"', '"HGDP00731"', '"HGDP00693"', '"HGDP00698"', '"HGDP00744"',
                     '"HGDP00732"', '"HGDP00679"', '"HGDP00686"', '"HGDP00675"', '"HGDP00726"', '"HGDP00733"',
                     '"HGDP00738"', '"HGDP00740"', '"HGDP00745"', '"HGDP00694"', '"HGDP00699"', '"HGDP00682"',
                     '"HGDP00687"')
Plestinian_Israel <- search_mongo(IDs_palestinian) %>%
  limpieza_total()

usethis::use_data(Plestinian_Israel, overwrite = TRUE)


# Papuan in Papua New Guinea HGDP

IDs_papuan <- c('"HGDP00549"', '"HGDP00551"', '"HGDP00556"', '"HGDP00553"', '"HGDP00550"', '"HGDP00555"')
Papuan_PapuaNewGuinea <- search_mongo(IDs_papuan) %>%
  limpieza_total()

usethis::use_data(Papuan_PapuaNewGuinea, overwrite = TRUE)


# Popuan Sepik in New Guinea HGDP

IDs_popuansepik <- c('"HGDP00544"', '"HGDP00540"', '"HGDP00545"', '"HGDP00541"', '"HGDP00546"',
                     '"HGDP00542"', '"HGDP00547"', '"HGDP00543"')
Popuan_Sepik_NewGuinea <- search_mongo(IDs_popuansepik) %>%
  limpieza_total()

usethis::use_data(Popuan_Sepik_NewGuinea, overwrite = TRUE)


# Pathan in Pakistan
IDs_pathan <- c('"HGDP00222"', '"HGDP00234"', '"HGDP00239"', '"HGDP00241"', '"HGDP00258"', '"HGDP00254"',
                '"HGDP00259"', '"HGDP00228"', '"HGDP00230"', '"HGDP00247"', '"HGDP00262"', '"HGDP00224"',
                '"HGDP00243"', '"HGDP00248"', '"HGDP00213"', '"HGDP00218"', '"HGDP00237"', '"HGDP00244"',
                '"HGDP00251"', '"HGDP00214"', '"HGDP00264"', '"HGDP00226"')
Pathan_Pakistan <- search_mongo(IDs_pathan) %>%
  limpieza_total()

usethis::use_data(Pathan_Pakistan, overwrite = TRUE)


# Pima in Mexico HGDP
IDs_Pima_SGCDP <- c('"HGDP01059"', '"HGDP01043"', '"HGDP01050"', '"HGDP01055"', '"HGDP01056"',
                    '"HGDP01037"', '"HGDP01051"', '"HGDP01057"', '"HGDP01041"', '"HGDP01053"',
                    '"HGDP01058"', '"HGDP01060"')
Pima_Mexico <- search_mongo(IDs_Pima_SGCDP) %>%
  limpieza_total()

usethis::use_data(Pima_Mexico, overwrite = TRUE)


# Russian in Russia HGDP

IDs_russian <- c('"HGDP00880"', '"HGDP00885"', '"HGDP00892"', '"HGDP00897"', '"HGDP00900"', '"HGDP00879"',
                 '"HGDP00881"', '"HGDP00886"', '"HGDP00893"', '"HGDP00898"', '"HGDP00901"', '"HGDP00882"',
                 '"HGDP00894"', '"HGDP00899"', '"HGDP00902"', '"HGDP00883"', '"HGDP00888"', '"HGDP00890"',
                 '"HGDP00895"', '"HGDP00889"', '"HGDP00891"', '"HGDP00896"', '"HGDP00884"')
Russian_Russia <- search_mongo(IDs_russian) %>%
  limpieza_total()

usethis::use_data(Russian_Russia, overwrite = TRUE)


# San in Namibia HGDP

IDs_san <- c('"HGDP01029"', '"HGDP00992"')
San_Namibia <- search_mongo(IDs_san) %>%
  limpieza_total()

usethis::use_data(San_Namibia, overwrite = TRUE)


# Sardinian in Italy HGDP

IDs_sardinian <- c('"HGDP00669"', '"HGDP00671"', '"HGDP01066"', '"HGDP01073"', '"HGDP00672"', '"HGDP01074"',
                   '"HGDP01062"', '"HGDP01067"', '"HGDP00666"', '"HGDP00673"', '"HGDP01070"', '"HGDP01075"',
                   '"HGDP01063"', '"HGDP01068"', '"HGDP00667"', '"HGDP00674"', '"HGDP01069"', '"HGDP01071"',
                   '"HGDP01064"', '"HGDP00668"', '"HGDP00670"', '"HGDP01065"', '"HGDP01072"', '"HGDP01077"')
Sardinian_Italy <- search_mongo(IDs_sardinian) %>%
  limpieza_total()

usethis::use_data(Sardinian_Italy, overwrite = TRUE)


# She in China HGDP

IDs_she <- c('"HGDP01327"', '"HGDP01334"', '"HGDP01328"', '"HGDP01330"', '"HGDP01329"', '"HGDP01331"',
             '"HGDP01336"', '"HGDP01332"')
She_China <- search_mongo(IDs_she) %>%
  limpieza_total()

usethis::use_data(She_China, overwrite = TRUE)


# Sindhi in Pakistan HGDP

IDs_sindhi <- c('"HGDP00183"', '"HGDP00210"', '"HGDP00169"', '"HGDP00171"', '"HGDP00165"', '"HGDP00177"',
                '"HGDP00189"', '"HGDP00191"', '"HGDP00185"', '"HGDP00192"', '"HGDP00197"', '"HGDP00205"',
                '"HGDP00173"', '"HGDP00181"', '"HGDP00201"', '"HGDP00206"', '"HGDP00167"', '"HGDP00179"',
                '"HGDP00187"', '"HGDP00199"', '"HGDP00163"', '"HGDP00175"')
Sindhi_Pakistan <- search_mongo(IDs_sindhi) %>%
  limpieza_total()

usethis::use_data(Sindhi_Pakistan, overwrite = TRUE)


# Surui in Brazil HGDP

IDs_surui <- c('"HGDP00843"', '"HGDP00832"', '"HGDP00837"', '"HGDP00849"', '"HGDP00838"', '"HGDP00845"')
Surui_Brazil <- search_mongo(IDs_surui) %>%
  limpieza_total()

usethis::use_data(Surui_Brazil, overwrite = TRUE)


# Tu in China HGDP

IDs_tu <- c('"HGDP01352"', '"HGDP01353"', '"HGDP01347"', '"HGDP01354"', '"HGDP01348"', '"HGDP01349"',
            '"HGDP01351"', '"HGDP01356"')
Tu_China <- search_mongo(IDs_tu) %>%
  limpieza_total()

usethis::use_data(Tu_China, overwrite = TRUE)


# Tujia in China HGDP

IDs_tujia <- c('"HGDP01097"', '"HGDP01100"', '"HGDP01101"', '"HGDP01099"', '"HGDP01102"', '"HGDP01103"',
               '"HGDP01096"', '"HGDP01104"')
Tujia_China <- search_mongo(IDs_tujia) %>%
  limpieza_total()

usethis::use_data(Tujia_China, overwrite = TRUE)


# Tuscan in Italy HGDP

IDs_tuscan <- c('"HGDP01162"', '"HGDP01167"', '"HGDP01164"', '"HGDP01169"', '"HGDP01166"', '"HGDP01161"')
Tuscan_Italy <- search_mongo(IDs_tuscan) %>%
  limpieza_total()

usethis::use_data(Tuscan_Italy, overwrite = TRUE)


# Uygur in China HGDP

IDs_uygur <- c('"HGDP01299"', '"HGDP01302"', '"HGDP01303"', '"HGDP01304"', '"HGDP01300"', '"HGDP01305"',
               '"HGDP01298"', '"HGDP01301"')
Uygur_China <- search_mongo(IDs_uygur) %>%
  limpieza_total()
Uygur_China <- readRDS("Uygur_China.rds")
Uygur_China <- NULL


# Xibo in China HGDP

IDs_xibo <- c('"HGDP01249"', '"HGDP01251"', '"HGDP01244"', '"HGDP01245"', '"HGDP01247"', '"HGDP01243"',
              '"HGDP01248"')
Xibo_China <- search_mongo(IDs_xibo) %>%
  limpieza_total()

usethis::use_data(Xibo_China, overwrite = TRUE)


# Yakut in Siberia HGDP

IDs_yakut <- c('"HGDP00948"', '"HGDP00950"', '"HGDP00955"', '"HGDP00962"', '"HGDP00967"', '"HGDP00963"',
               '"HGDP00968"', '"HGDP00949"', '"HGDP00964"', '"HGDP00969"', '"HGDP00945"', '"HGDP00952"',
               '"HGDP00957"', '"HGDP00960"', '"HGDP00965"', '"HGDP00946"', '"HGDP00953"', '"HGDP00958"',
               '"HGDP00961"', '"HGDP00966"', '"HGDP00947"', '"HGDP00954"', '"HGDP00959"')
Yakut_Siberia <- search_mongo(IDs_yakut) %>%
  limpieza_total()

usethis::use_data(Yakut_Siberia, overwrite = TRUE)


# Yi in China

IDs_yi <- c('"HGDP01181"', '"HGDP01186"', '"HGDP01187"', '"HGDP01182"', '"HGDP01183"', '"HGDP01184"', '"HGDP01180"',
            '"HGDP01185"')
Yi_China <- search_mongo(IDs_yi) %>%
  limpieza_total()

usethis::use_data(Yi_China, overwrite = TRUE)


# Yoruba in Nigeria HGDP

IDs_yoruba <- c('"HGDP00943"', '"HGDP00924"', '"HGDP00929"', '"HGDP00931"', '"HGDP00920"', '"HGDP00944"',
                '"HGDP00925"', '"HGDP00937"', '"HGDP00926"', '"HGDP00933"', '"HGDP00938"', '"HGDP00940"',
                '"HGDP00934"', '"HGDP00939"', '"HGDP00941"', '"HGDP00930"', '"HGDP00935"', '"HGDP00942"')
Yoruba_Nigeria <- search_mongo(IDs_yoruba) %>%
  limpieza_total()

usethis::use_data(Yi_China, overwrite = TRUE)






