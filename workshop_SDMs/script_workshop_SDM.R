####################################################################################################
# Workshop Species Distribution Modelling (SDM), STOWA CoP ecologische data-analyse                #
# Dit script dient als een eerste kennismaking met SDMs. SDMs zijn een handig instrument           #
# bij het voorspellen van de distributie van plant- en diersoorten. Deze modellen zijn bijvoorbeeld#
# populair in onderzoek dat zich focust op het effect van klimaatverandering op de verspreidings-  #
# areaal van een bepaalde soort. In dit script zal de ruimtelijke distributie van de               #
# alpenwatersalamander in Nederland, nu en over een aantal decennia, worden gemodelleerd.          #
####################################################################################################
# Eerst hebben we de benodigde software nodig, verpakt in zogeheten 'packages'. In deze packages
# bevinden zich de functies die we later nodig hebben. In deze workshop staat het biomod2
# pakket centraal. Dit pakket biedt de opties waarmee we kunnen modelleren.
# install.packages('biomod2') # Om het pakket de eerste keer te kunnen gebruiken moet het eerst
# geinstalleerd worden.
library(biomod2) # het pakket waarmee we kunnen modelleren
# install.packages('tidyverse') # wederom, bij eerste gebruik, eerst installeren
library(tidyverse)
# install.packages('terra') # bij eerste gebruik eerst installeren
library(terra)
# install.packages('tidyterra') # bij eerste gebruik eerst installeren
library(tidyterra)
# install.packages('raster') # bij eerste gebruik eerst installeren
library(raster)
# install.packages('dismo') # bij eerste gebruik eerst installeren
library(dismo)
# install.packages('ggplot2') # bij eerste gebruik eerst installeren
library(ggplot2)

# En, belangrijk: we moeten onze 'working directory' definiëren, dit is het pad op onze pc
# waar we de benodigde bestanden kunnen ophalen. Dit kan met de 'setwd()' functie. Let op:
# R kan dit alleen 'lezen' als je '/' gebruikt en niet '\'. Het makkelijkst in dit geval is het
# opslaan van de benodigde data onder een makkelijk te vinden mapje, zoals bijvoorbeeld
# 'documenten', en hierbinnen dan een nieuwe folder te creëren met de naam "SDM".
# Hierbinnen is het dan handig om een sub-folder te hebben met de naam "input_data". Hier
# kan de gedownloade data in worden gezet, zodat deze zodadelijk makkkelijk kan worden ingeladen.
# Working directory instellen:
setwd("C:/Documents/SDM") # pad naar SDM map, pad naar map met data hier kopiëren
####################################################################################################
# Voor het opstellen van SDMs hebben we twee 'basisingrediënten' nodig:                            #
# - Waarnemingsdata: locatiegebonden punten van de soort waarvoor we een model willen opstellen    #
# - Omgevingsvariabelen: we willen graag inzicht krijgen in de relatie tussen het wel of niet      #
#   voorkomen van een soort en een selectie van omgevingsvariabelen. Dit kunnen zowel              #
#   klimatologische variabelen (als temperatuur, neerslag) zijn, maar het is ook mogelijk          #
#   om hier categorische variabelen zoals landgebruikstype in te stoppen.                          #
####################################################################################################
# Allereerst de waarnemingsdata. In dit script zullen we werken met  Nederlandse 
# waarnemingsdata voor alpenwatersalamander in de afgelopen 15 jaar.
# inladen en data opslaan onder de naam 'waarnemingsdata':
waarnemingsdata <- read.csv("./input_data/data_alpenwatersalamander.csv")
######################################################################################################
# Naast de waarnemingsdata hebben we ook nog omgevingsvariabelen nodig. In dit geval zullen
# we het houden bij klimatologische variabelen. We zullen werken met gemiddelde
# jaarlijkse temperatuur, seizoensgebonden variatie in temperatuur, jaarlijkse hoeveelheid
# en seizoensgebonden variatie in neerslag. Deze vier variabelen beschrijven de belangrijkste
# klimatologische kenmerken van een locatie. Hiernaast laden we ook de gemiddelde temperatuur
# van het koudste kwartaal in, dit is ook een belangrijk klimatologisch kenmerk voor amfibieën.
# Inladen van deze data voor Nederland:
omgevingsdata <- stack(c(gem_T = raster("./input_data/gem_temp.tif"), 
                         var_T = raster("./input_data/var_temp.tif"),
                         min_kwart_T = raster("./input_data/koudste_kwartaal.tif"),
                         gem_N = raster("./input_data/gem_neerslag.tif"),
                         var_N = raster("./input_data/var_neerslag.tif")))
omgevingsdata # overzichtje van kenmerken omgevingsdata.
plot(omgevingsdata) # visualisatie van de input-omgevingsdata.
#####################################################################################################
# Nu hebben we de belangrijkste input-data klaarstaan. We kunnen nu beginnen met het bouwen
# van modellen. Hiervoor gaan we eerst de data opslaan in een format waar de software mee kan 
# rekenen. Een belangrijk punt in deze stap is het genereren van zogeheten pseudo-absenties (PAs).
# SDMs hebben namelijk ook 'niet-waarnemingen' nodig. In dit geval hebben we alleen aanwezigheids-
# data, dus we zullen deze 'achtergrondpunten' zelf moeten genereren. We zullen dit nu doen
# met basis-instellingen, maar biomod2 biedt een hoop opties om je set PAs te tweaken.
#####################################################################################################
modelData <- BIOMOD_FormatingData(resp.var = waarnemingsdata[,2], # aanwezigheidsdata
                                  resp.xy = waarnemingsdata[,c("lengtegraad", "breedtegraad")], # coördinaten van waarnemingen
                                  expl.var = omgevingsdata,   # omgevingsvariabelen
                                  resp.name = 'alpenwatersalamander', # naam van output-map
                                  PA.strategy = "random",     # standaard strategie voor opstellen PAs
                                  PA.nb.rep = 1,              # aantal PA sets, standaard is 1
                                  PA.nb.absences = 13000,     # aantal PAs in de set
                                  filter.raster = TRUE)       # waarnemingsdata wordt gefilterd
                                                              # bij meerdere records in 1 cel.
######################################################################################################
# Met deze data kunnen we nu daadwerkelijk de SDMs op gaan stellen!
# In deze stap zullen we de geformatteerde data uit de vorige stap gebruiken. Hiernaast
# zullen we de data in deze stap splitten in vier blokken. Hiervan worden drie blokken gebruikt
# voor het calibreren van het model, terwijl het laatste blok wordt gebruikt om de voorspelling
# van het model tegen te kunnen evalueren.
modOptions <- BIOMOD_ModelingOptions()
SDMs <- BIOMOD_Modeling(bm.format = modelData, # hierboven geprepareerde data
                        bm.options = modOptions,
                        modeling.id = "alpenwatersalamander_huidig", # 'tag' voor deze modelrun
                        models = c('GLM', 'CTA'), # te gebruiken algoritme(s)
                        CV.strategy = "kfold", # strategie voor 'kruisvalidatie'
                        CV.nb.rep = 2,         # aantal herhalingen kruisvalidatie
                        CV.k = 4,              # aantal blokken waarin we data splitten
                        metric.eval = c("TSS", "ROC"), # parameters waarmee modelprestatie wordt beoordeeld
                        var.import = 2) # aantal permutaties omgevingsvariabelen om belang 
                                        # variabele in te kunnen schatten.
######################################################################################################
# Met deze uitkomst kunnen we nu de prestatie van het model gaan evalueren en de resultaten
# visualiseren. Allereerst de evaluatie op basis van TSS en ROC. Hierbij kan TSS ('True
# Skill Statistic') een waarde aannemen tussen -1 en 1, waarbij 1 een perfecte voorspelling
# van de werkelijkheid betekent, 0 aangeeft dat het model even goed is als een willekeurige
# voorspelling, en -1 een totale discrepantie tussen werkelijkheid en voorspellingen weergeeft.
# ROC ('Receiver Operating Characteristic') kan tussen 0 en 1 liggen. Hierbij is wederom 1
# een model dat gelijk is aan de werkelijkheid terwijl 0.5 een prestatie gelijk aan willekeurig
# betekent, en 0 betekent dat je model het elke keer verkeerd heeft.
summary(get_evaluations(SDMs))
# deze data kunnen we ook visualiseren in een boxplot.
bm_PlotEvalMean(bm.out = SDMs,     # onze SDMs gaan erin
                group.by = 'algo') # ze worden gescheiden op algoritme
# In dit plaatje zie je ook duidelijk de verschillen tussen prestaties van de 
# verschillende algorithmen (GLM en CTA).
# We kunnen op een vegelijkbare manier ook de invloed van de verschillende variabelen op 
# de voorspellingen van de modellen analyseren.
get_variables_importance(SDMs)
# Ook deze data kan worden gevisualiseerd: 
bm_PlotVarImpBoxplot(bm.out = SDMs,                    # onze modellen erin
                     group.by = c('algo', 'expl.var', 'PA')) # gesplit of basis van variabele en algoritme
# Deze plots laten zien hoeveel van de variatie in de modelberekeningen verklaard kunnen worden
# door variatie in de omgevingsvariabele in kwestie.

# Als laatste stap kan ook de geschiktheid van habitat als functie van waarden van de verschillende
# variabelen worden gevisualiseerd in zogeheten 'response curves'.
rc <- bm_PlotResponseCurves(bm.out = SDMs,   # onze modellen
                      models.chosen = 'all', # alle berekende modellen worden meegenomen
                      fixed.var = 'mean')    # variabelen die niet worden veranderd worden gemiddeld houden
# dit plaatje is een beetje moeilijk te interpreteren. Hieronder gaan we een overzichtelijker
# plot maken
rc_data <- rc$tab # data uit de respons-curves halen
# we willen kunnen filteren op algoritme:
rc_data$alg <- substr(rc_data$pred.name, start = nchar(as.character(rc_data$pred.name))-2, stop = 
                        nchar(as.character(rc_data$pred.name))) 
# en op run:
rc_data$run <- substr(rc_data$pred.name, start = nchar(as.character(rc_data$pred.name))-7, stop = 
                        nchar(as.character(rc_data$pred.name))-4)
# Duidelijker plot maken:
ggplot(data = rc_data, aes(x = expl.val, y = pred.val, color = run)) +
  theme_bw() +
  geom_line() +
  facet_grid(alg ~ expl.name, scales = "free_x")
# In deze figuur zie je op de y-as de geschiktheid van het leefgebied in relatie tot 
# de waarden van de omgevingsvariabelen, welke uitgesmeerd zijn over de x-as. Wederom
# duideljke verschillen in de relaties (en de vorm van de relatie) die de verschillende
# algoritmen voorspellen.
####################################################################################################
# We kunnen de modellen ook projecteren en de output visualiseren:
projectie_Nederland <- BIOMOD_Projection(bm.mod = SDMs, # onze SDMs
                                         proj.name = "Nederland",
                                         new.env = omgevingsdata,# klimaatdata voor NL
                                         selected.models = 'all', # voor alle modellen
                                         metric.binary = 'all', # binariseren met zelfde parameters als modelevaluatie
                                         metric.filter = 'all')# filteren met zelfde parameters als modelevaluatie
# deze projectie kan gevisualiseerd worden. Vergroten kan met 'zoom' knopje in plot paneel.
plot(projectie_Nederland, plot.output = 'facet') # Verschillende projecties naast elkaar,
# Hierbij geven de lichtere kleuren een hogere geschiktheid voor onze soort weer.
# Dit kan ook binair worden weergegeven (0 = ongeschikt, 1 = geschikt):
ranges <- get_predictions(projectie_Nederland,  # eerder opgestelde projecties
                         metric.binary = 'TSS',# binariseren op basis van TSS
                         model.as.col = TRUE)  
plot(ranges, main = "geschiktheid habitat",
     plot.output = 'facet')
# Het resultaat is een kaart met in het groen gebieden welke door het model aan worden gestipt
# als mogelijk geschikt voor alpenwatersalamander en in het grijs gebieden welke niet geschikt
# worden geacht. Let wederom op de verschillen tussen de rekenmethoden.
# Dit waren de basics van het biomod2 pakket en het opstellen van SDMs. Bij behoefte aan meer 
# verdieping, en tijd over, zijn er hieronder nog een paar optionele stappen te vinden.
#############################################################################################################
# Een populaire volgende stap is het opstellen van een ensemble model. Dit model combineert
# de verschillende algoritmen in een overkoepelend model, waarmee er minder nadruk komt 
# te liggen op de zwakten van de verschillende rekenmethoden, en de sterke punten gebundeld
# kunnen worden.
ensemble_SDM <- BIOMOD_EnsembleModeling(bm.mod = SDMs, # onze individuele modellen
                                        models.chosen = 'all', # alle modellen meenemen
                                        em.by = 'all', # combineren tot één overkoepelend model
                                        em.algo = 'EMwmean', # beter presterende modellen krijgen meer gewicht
                                        metric.select = 'ROC', # op basis van ROC worden modellen wel of niet meegenomen
                                        metric.select.thresh = 0.75, # ROC moet minimaal 0.75 zijn voor opname ensemble model
                                        metric.eval = c('TSS', 'ROC'), # evaluatie op basis van ROC en TSS
                                        EMwmean.decay = 'proportional', # minder presterende modellen krijgen propoortioneel minder gewicht
                                        var.import = 1, # aantal permutaties omgevingsvariabelen om belang 
                                        # variabele in te kunnen schatten.
                                        do.progress = TRUE) # balkje dat laat zien hoever R is.
###################################################################################################################################
# Dit ensemble model kan weer op dezelfde manier worden geanalyseerd worden als de individuele modellen".
# Wederom gebeurt dit aan de hand van ROC en TSS.
get_evaluations(ensemble_SDM)
# deze data kunnen we ook visualiseren in een boxplot, al is het maar één datapunt dus 
# het ziet er een beetje vreemd uit.
bm_PlotEvalMean(bm.out = ensemble_SDM,     # ons ensemble model gaat erin
                group.by = 'algo')
# ook het belang van de variabelen in het ensemble model kan geanalyseerd worden:
get_variables_importance(ensemble_SDM)
# Ook deze data kan worden gevisualiseerd: 
bm_PlotVarImpBoxplot(bm.out = ensemble_SDM, # ons ensemble model er weer in
                     group.by = c('algo', 'expl.var', 'merged.by.PA')) 
# En ook voor het ensemble model kunnen weer respons-curves worden opgesteld, via een
# vergelijkbare werkwijze als voorheen:
bm_PlotResponseCurves(bm.out = ensemble_SDM, # ons ensemble-model
                      models.chosen = 'all', # alle berekende modellen worden meegenomen
                      fixed.var = 'mean')    # niet getweakte variabelen weer op gemiddelde waarde
# In deze figuur zie je wederom op de y-as de geschiktheid van het leefgebied in relatie tot 
# de waarden van de omgevingsvariabelen, welke uitgesmeerd zijn over de x-as.
# Ook van dit ensemble model kunnen we weer een projectie maken op continue en binaire schaal:
# eerst de projectie:
ensemble_projectie <- BIOMOD_EnsembleForecasting(bm.em = ensemble_SDM, # ons model
                                                 proj.name = 'ensemble_NL', # naam projectie
                                                 new.env = omgevingsdata, # huidige data
                                                 models.chosen = 'all', 
                                                 metric.binary = 'all',
                                                 metric.filter = 'all')
# deze projectie kan wederom gevisualiseerd worden:
plot(ensemble_projectie) 
# Hierbij geven de lichtere kleuren een hogere geschiktheid voor onze soort weer.
# Ook deze projectie kan binair worden weergegeven (0 = ongeschikt, 1 = geschikt):
ensemble_range <- get_predictions(ensemble_projectie,  # ensemble projectie
                                  metric.binary = 'TSS', # binariseren op basis van TSS
                                  model.as.col = TRUE)  
plot(ensemble_range, main = "geschiktheid habitat")
# Het resultaat is weer een kaart met in het groen gebieden welke door het model aan worden 
# gestipt als mogelijk geschikt voor alpenwatersalamander en in het grijs gebieden welke
# niet geschikt worden geacht. In theorie zou dit plaatje overeen moeten komen met de huidige
# verspreiding van alpenwatersalamander in Nederland.
###############################################################################################
# Nog een populaire toepassing van (ensemble) SDMs is het doorrekenen van klimaatscenario's
# en het effect hiervan op de distributie van een soort beoordelen. Dit is het laatste 
# punt dat uitgewerkt wordt in dit script.
# Inladen projectie voor 2060:
omgevingsdata_2060 <- stack(c(gem_T = raster("./input_data/gem_temp_NL_toekomst.tif"), 
                            var_T = raster("./input_data/var_temp_NL_toekomst.tif"),
                            min_kwart_T = raster("./input_data/koudste_kwartaal_NL_toekomst.tif"),
                            gem_N = raster("./input_data/gem_neerslag_NL_toekomst.tif"),
                            var_N = raster("./input_data/var_neerslag_NL_toekomst.tif")))
omgevingsdata_2060 # kenmerken van deze data
plot(omgevingsdata_2060) # en gevisualiseerd
###################################################################################################
# Nu kunnen we de modellen, in dit geval het ensemble model, projecteren naar deze omstandigheden
# en visualiseren wat er volgens dit model verandert als het toekomstige klimaat het 
# traject van deze projectie volgt.
ensemble_projectie_2060 <- BIOMOD_EnsembleForecasting(bm.em = ensemble_SDM, # ons model
                                                      proj.name = 'ensemble_NL_2060', # naam projectie
                                                      new.env = omgevingsdata_2060, # data 2080
                                                      models.chosen = 'all', 
                                                      metric.binary = 'all',
                                                      metric.filter = 'all')

# deze projectie kunnen we weer visualiseren op de kaart, zowel op continue als binaire schaal:
plot(ensemble_projectie_2060) 
# Hierbij geven de lichtere kleuren een hogere geschiktheid voor onze soort weer.
# Ook deze projectie kan binair worden weergegeven (0 = ongeschikt, 1 = geschikt):
ensemble_range_2060 <- get_predictions(ensemble_projectie_2060,  # ensemble projectie
                                       metric.binary = 'TSS', # binariseren op basis van TSS
                                       model.as.col = TRUE)  
plot(ensemble_range_2060, main = "geschiktheid habitat 2060")
# Het resultaat is weer een kaart met in het groen gebieden welke door het model aan worden 
# gestipt als mogelijk geschikt voor alpenwatersalamander en in het grijs gebieden welke
# niet geschikt worden geacht.
# Het verschil in range kan ook mooi worden gevisualiseerd met dit pakket:
range_verschil <- BIOMOD_RangeSize(ensemble_range, ensemble_range_2060)
range_verschil$Compt.By.Models # overzicht hoeveel pixels verloren zijn gegaan, hoeveel er
                               # stabiel ongeschikt en geschikt zijn gebleven en hoeveel
                               # er bij zijn gekomen.
# dit kan ook weer worden gevisualiseerd:
plot(range_verschil$Diff.By.Pixel)
# waarbij '-2' verlies van geschikt habitat laat zien, '0' stabiel gebleven (on)geschikt habitat
# en '1' nieuw geschikt habitat laat zien.
###############################################################################################
# Dit waren in een notendop de algemene toepassingen van Species Distribution Models.