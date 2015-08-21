data(Owls, package="glmmADMB")
Owls <- transform(Owls, ArrivalTime=c(scale(ArrivalTime,scale=FALSE)),
                  NCalls=SiblingNegotiation)
tmbz <- glmmTMB(NCalls ~ FoodTreatment + ArrivalTime + (ArrivalTime|Nest) + (1|SexParent),
                ziformula = ~1 + (1|Nest) + (1|SexParent),
                data=Owls, family=poisson(link="log"))
object <- tmbz

--------------------------------------------------------------------------------

BEFORE

> ranef(object)

$`conditional model`
     AutavauxTV      AutavauxTV          Bochet          Bochet     Champmartin     Champmartin
  -1.799535e-01    1.483289e-01    7.185302e-02   -5.801783e-03   -7.486459e-02   -2.579779e-02
        ChEsard         ChEsard        Chevroux        Chevroux CorcellesFavres CorcellesFavres
   3.062864e-01   -1.143882e-01    1.703923e-01   -6.550542e-02   -2.394259e-01    1.511229e-01
       Etrabloz        Etrabloz           Forel           Forel          Franex          Franex
  -2.357961e-01    1.018474e-01   -3.148762e-01    1.066159e-01   -5.661592e-02    5.782232e-04
           GDLV            GDLV      Gletterens      Gletterens         Henniez         Henniez
  -4.772949e-02    5.345830e-02    4.275284e-02   -1.224408e-01   -1.400744e-01    5.948369e-02
          Jeuss           Jeuss     LesPlanches     LesPlanches          Lucens          Lucens
  -2.290793e-01   -3.311348e-02    1.154527e-02    3.930111e-02    7.030649e-02   -7.046802e-02
          Lully           Lully         Marnand         Marnand          Montet          Montet
   1.020255e-01   -9.184923e-02   -5.371682e-03   -4.662523e-02    2.384561e-01   -5.383558e-02
         Murist          Murist          Oleyes          Oleyes         Payerne         Payerne
   2.062517e-01   -1.067604e-01    1.817386e-01   -4.262321e-02    3.772022e-02    5.947097e-02
         Rueyes          Rueyes           Seiry           Seiry           Sevaz           Sevaz
   4.533439e-02   -1.005719e-01    1.534573e-01    8.676384e-02   -2.610783e-01    2.725364e-02
        StAubin         StAubin            Trey            Trey        Yvonnand        Yvonnand
  -3.651268e-01    1.391811e-01    3.175307e-01   -6.989095e-02    3.339912e-01   -2.898687e-02
         Female            Male
   4.798218e-10    2.403581e-10

$zero_inflation
     AutavauxTV          Bochet     Champmartin         ChEsard        Chevroux CorcellesFavres
     0.16723405      0.82768972      0.79057668     -0.88584882      1.02979235     -0.94102795
       Etrabloz           Forel          Franex            GDLV      Gletterens         Henniez
     0.27518184      0.92063539     -0.04222336     -0.46151325      0.66418950      0.98158951
          Jeuss     LesPlanches          Lucens           Lully         Marnand          Montet
     0.75281327      1.18174261     -0.04120996     -0.51550414     -0.09657730     -0.60817225
         Murist          Oleyes         Payerne          Rueyes           Seiry           Sevaz
     0.20597202     -0.09263899      0.14896882     -0.24771061     -0.85286112     -0.47071429
        StAubin            Trey        Yvonnand          Female            Male
    -0.06259358     -0.93979375     -1.07230524      0.16377871     -0.12625047

--------------------------------------------------------------------------------

m <- matrix(ranef(tmbz)[[1]][1:54], ncol=2, byrow=TRUE)
gprVar <- object$modelInfo$grpVar
rownames(m) <-

AFTER

> ranef(object)

$`conditional model`

