# Lab29_02_24
Modificació del codi d'un pou quadrat per a obtenir un pou triangular.

Per a analitzar el cas en que el potencial no té un efecte significant en la zona de la barrera, cal definir el potencial de la següent manera:
def getV(x):
    if (abs(x)<(Lpou/2)):
       potvalue = 0.0 +Velectr(x)
    else:
       potvalue = 4.0 
    return potvalue
    
En cas de voler analitzar el cas en que el potencial presenta un efecte en la zona de la barrera, cal definir el potencial de la següent manera:
def getV(x):
    if (abs(x)<(Lpou/2)):
       potvalue = 0.0 +Velectr(x)
    else:
       potvalue = 4.0 +Velectr(x)
    return potvalue
