from psychopy import visual, core, event, gui
import glob
import random
import pandas as pd
#-------------

#Define pop-up window at the start of the experiment
box = gui.Dlg(title = "Smerte og hukommelse")
box.addField("ID*:")
box.addField("Navn*:")
box.addField("Alder")
box.addField("Køn:", choices = ["Vælg","Mand", "Kvinde"])
box.addField("Hvor mange kopper kaffe har du indtaget i dag:", choices = ["Vælg", "Intet kaffe", "Én kop", "To kopper", "Tre kopper", "Mere end tre kopper"])
box.addField("Har du brugt smertestillende midler, indenfor de seneste 24 timer?:", choices = ["Vælg", "Ja", "Nej"])

#Show pop-up box
box.show()

# Get the data
if box.OK:
    ID = box.data[0]
    Navn = box.data[1]
    Alder = box.data[2]
    Køn = box.data[3]
    Kaffe = box.data[4]
    Smertestillende = box.data[5]
elif box.Cancel:
    core.quit()

#Define window
win = visual.Window(fullscr = True, color = "white", units="pix")

#Define logfile and the colums within
logfile = pd.DataFrame(columns = ["ID", "Name", "Age", "Gender", "Coffee", "Painkillers", "Text", "Stimulus", "Question", "Perception", "Objective"])

#Samtykke
samtykke = visual.ImageStim(win, image = "Python materials/Samtykke.PNG", size=(1680,1050))
samtykke.draw()
win.flip()
key = event.waitKeys(keyList = ["space", "q"])
if key[0] == "q":
    core.quit()

#Introduction 1: "Velkommen"
velkommen = visual.ImageStim(win, image = "Python materials/Velkommen.PNG", size=(1680,1050))
velkommen.draw()
win.flip()
key = event.waitKeys(keyList = ["space", "q"])
if key[0] == "q":
    core.quit()

#Introduction 2: "Stød-skala"
stødskala = visual.ImageStim(win, image = "Python materials/Stødskala.PNG", size=(1680,1050))
stødskala.draw()
win.flip()
key = event.waitKeys(keyList = ["return", "q"])
if key[0] == "q":
    core.quit()

#Introduction 2: "Klar"
klar = visual.ImageStim(win, image = "Python materials/Klar.PNG", size=(1680,1050))
klar.draw()
win.flip()
key = event.waitKeys(keyList = ["space", "q"])
if key[0] == "q":
    core.quit()


#Define stories
B1 = visual.TextStim(win,'''Det var en varm og dejlig solskinsdag, jeg var lige kommet hjem og havde krammet min dejlige kæreste, da jeg fik det mest fortryllende kys nogensinde.''', height = 25, color = "black")
B2 = visual.TextStim(win,"Jeg flirter lidt med den dejlige mand jeg har mødt. Jeg kan simpelthen ikke lade være, han er helt og aldeles himmelsk, og så har han den sødeste og mest trofaste hund.", height = 25, color = "black")
B3 = visual.TextStim(win,"Et af de mest opfriskende minder, som gang på gang gør mig lykkelig. Er af min søde og dejlige lille niece, der gynger ude i haven på en varm solskinsdag.", height = 25, color = "black")
N1 = visual.TextStim(win,"Jeg skulle bo på et klassisk hotel, mens mit almindelige hus blev renoveret. Hotellet var rødt og blev bygget i 1999, og har udsigt ud over den enorme by.", height = 25, color = "black")
N2 = visual.TextStim(win,"Samtalen hos den foreløbige chef, havde været lang og vi snakkede kun om officielle emner, som ikke påvirkede borgerne. De lange samtaler var ved at blive dagligdag for mig.", height = 25, color = "black")
N3 = visual.TextStim(win,"Jeg havde altid været lav og tynd i hele gymnasietiden, men efter at jeg havde været ude på den lange internationale rejse, var der sket en klar ændring.", height = 25, color = "black")
S1 = visual.TextStim(win,"Jeg kom ind på hospitalet, det var forfærdeligt. Min arm var blevet skoldet og smerten blev ved med at føltes varm og sviende, det havde været en grusom oplevelse.", height = 25, color = "black")
S2 = visual.TextStim(win,"Smerten føltes jagende og sviende, da jeg vågnede op midt i operationen. De var igang med at skære mit venstre bryst op. Jeg frygtede for livet - tænkte at de ville dræbe mig.", height = 25, color = "black")
S3 = visual.TextStim(win,"Mit hjerte bankede og jeg havde en frygtelig tom fornemmelse i maven efter den grufulde og ødelæggende oplevelse. Det var forfærdeligt og jeg følte mig kvalt i situationen.", height = 25, color = "black")

#pause
pause = visual.TextStim(win,"Du vil nu blive præsenteret for nogle sætninger", height = 25, color = "black", bold=True)


#Define spørgsmål
QB1 = visual.TextStim(win,''' Hvad skete der da personen var kommet hjem?


1. Han krammede sin dejlige kæreste og fik et fortryllende kys.

2. Han mødte sin venlige kæreste der krammede ham

3. Han kyssede sin kæreste og de gik ud i det fine vejr ''', height = 20, color = "black")
Q2B1 = visual.TextStim(win,''' Hvem mødte personen da han kom hjem?


1. Hans kæreste

2. Hans hund

3. Hans søster. ''', height = 20, color = "black")
QB2 = visual.TextStim(win,''' Hvorfor flirter personen med den dejlige mand?


1. Fordi han var sød og dejlig

2. Fordi han havde hjulpet hende igennem en svær tid

3. Fordi han var himmelsk og havde en sød og trofast hund. ''', height = 20, color = "black")
Q2B2 = visual.TextStim(win,''' Hvilket ord blevet nævnt i historien?

1. Dejlig

2. Solskin

3. Kvinde ''', height = 20, color = "black")
QB3 = visual.TextStim(win,''' Hvilket minde gjorde personen lykkelig gang på gang?


1. Mindet om en sød og dejlig niece, der gynger i solskinsvejr 

2. Mindet om en sød og dejlig niece der leger ude i solskinsvejret 

3. Mindet om den søde niece som hoppede i vandpytter ''', height = 20, color = "black")
Q2B3 = visual.TextStim(win,''' Hvad gjorde niecen i historien?


1. Gyngede

2. Kravlede

3. Rutsjede ''', height = 20, color = "black")

QN1 = visual.TextStim(win,''' Hvorfor boede personen på et hotel?


1. Det var et flot hotel fra 1999

2. Grundet renovering af eget hus

3. Han skulle se den gode udsigt''', height = 20, color = "black")
Q2N1 = visual.TextStim(win,''' Hvilken farve var hotellet?


1. Rødt

2. Blåt

3. Gult
 ''', height = 20, color = "black")
QN2 = visual.TextStim(win,''' Hvad handlede samtalen i sætningen om?


1. Hvem der skulle være den nye chef

2. Hvordan samtalerne var blevet dagligdag.

3. Emner der ikke påvirkede borgerne
 ''', height = 20, color = "black")
Q2N2 = visual.TextStim(win,''' Hvilket af følgende ord blev ikke nævnt i historien?


1. Chef

2. Almindelig

3. Kedelig ''', height = 20, color = "black")
QN3 = visual.TextStim(win,''' Hvad var der sket med personen efter at han havde været ude på sin rejse?


1. Han havde fået en masse internationale kontakter

2. Han var blevet verdenskendt og havde fået en masse venner

3. Han havde taget en del kilo på og blevet højere.
 ''', height = 20, color = "black")
Q2N3 = visual.TextStim(win,''' Hvordan så personen ud i gymnasietiden? 


1. Stor og flot

2. Tynd og lav

3. Høj og tynd
 ''', height = 20, color = "black")

QS1 = visual.TextStim(win,''' Hvad kunne personen i historien have været ude for?


1. Han havde fået tredje-grads forbrændinger

2. Hans far var lige død og skulle ind og se ham for sidste gang

3. Hans kæreste var blevet indlagt efter et slagtilfælde ''', height = 20, color = "black")
Q2S1 = visual.TextStim(win,''' Hvor foregår historien?


1. På et hospital 

2. I et hus

3. På skadestuen ''', height = 20, color = "black")
QS2 = visual.TextStim(win,''' Hvorfor tror du at personen skulle opereres?


1. Dårlige knæ

2. Problemer med hjertet

3. Fået en kniv i maven  ''', height = 20, color = "black")
Q2S2 = visual.TextStim(win,''' Hvad var de igang med da personen vågende op?


1. At skære i personens højre bryst

2. At skære i personens venstre bryst ''', height = 20, color = "black")
QS3 = visual.TextStim(win,''' Hvilken situation kunne passe på historien?


1. Personens bror faldt om på gaden

2. Personen faldt og slog sit knæ

3. Personen væltede på sin cykel. 
 ''', height = 20, color = "black")
Q2S3 = visual.TextStim(win,''' Hvad havde oplevelsen været?

1. Tom og ødelæggende

2. Grufuld og forfærdelig

3. Grufuld og ødelæggende ''', height = 20, color = "black")

#Define text function
def text(x):
    x.draw()
    win.flip()
    core.wait(15)

#Randomize stories
list1 = [[B1,"1", QB1, Q2B1, 1, 9],[B2,"3",QB2, Q2B2, 2, 2],[B3,"1",QB3, Q2B3, 3, 5],[N1,"2",QN1, Q2N1, 4, 7],[N2,"3",QN2, Q2N2, 5, 3],[N3,"3",QN3, Q2N3, 6, 5],[S1,"1",QS1, Q2S1, 7, 8],[S2,"2",QS2, Q2S2, 8, 4],[S3,"1",QS3, Q2S3, 9, 6]]
list2 = [[B1,"1", QB1, Q2B1, 1, 7],[B2,"3",QB2, Q2B2, 2, 3],[B3,"1",QB3, Q2B3, 3, 5],[N1,"2",QN1, Q2N1, 4, 8],[N2,"3",QN2, Q2N2, 5, 4],[N3,"3",QN3, Q2N3, 6, 6],[S1,"1",QS1, Q2S1, 7, 9],[S2,"2",QS2, Q2S2, 8, 2],[S3,"1",QS3, Q2S3, 9, 5]]
list3 = [[B1,"1", QB1, Q2B1, 1, 6],[B2,"3",QB2, Q2B2, 2, 8],[B3,"1",QB3, Q2B3, 3, 4],[N1,"2",QN1, Q2N1, 4, 5],[N2,"3",QN2, Q2N2, 5, 9],[N3,"3",QN3, Q2N3, 6, 2],[S1,"1",QS1, Q2S1, 7, 5],[S2,"2",QS2, Q2S2, 8, 7],[S3,"1",QS3, Q2S3, 9, 3]]
list4 = [[B1,"1", QB1, Q2B1, 1, 7],[B2,"3",QB2, Q2B2, 2, 9],[B3,"1",QB3, Q2B3, 3, 2],[N1,"2",QN1, Q2N1, 4, 5],[N2,"3",QN2, Q2N2, 5, 6],[N3,"3",QN3, Q2N3, 6, 2],[S1,"1",QS1, Q2S1, 7, 5],[S2,"2",QS2, Q2S2, 8, 8],[S3,"1",QS3, Q2S3, 9, 3]]

random.shuffle(list1)
random.shuffle(list2)
random.shuffle(list3)
random.shuffle(list4)

#Present text
for i in list4:
    pause.draw()
    win.flip()
    core.wait(3)
    text(i[0])
    i[2].draw()
    win.flip()
    key = event.waitKeys(keyList = ["1", "2","3","4","q"])
    if key[0] == i[1]:
        answer = "Correct"
    elif key[0] == "q":
        core.quit()
    elif key[0] != i[1]:
        answer = key[0]
    pic = visual.ImageStim(win, image = "Python materials/Scale.PNG", size=(1680,1050))
    pic.draw()
    win.flip()
    scalekey = event.waitKeys(keyList = ["dollar","1","2","3","4","5","6","7","8","9","0","q"])
    if scalekey[0] == "q":
        core.quit()
    elif scalekey[0] == "0":
        perception = "10"
    elif scalekey[0] == "dollar":
        perception = "0"
    elif scalekey[0] != "dollar":
        perception = scalekey[0]
    i[3].draw()
    win.flip()
    key = event.waitKeys(keyList = ["1", "2","3","4","q"])
    print(answer)
    if i[4] <= 3:
        Stimulus = "Comfortable"
    elif i[4] <= 6:
        Stimulus = "Neutral"
    else:
        Stimulus = "Painful"
    logfile = logfile.append({
            "ID": ID,
            "Name": Navn,
            "Age": Alder,
            "Gender": Køn,
            "Coffee": Kaffe,
            "Painkillers": Smertestillende,
            "Text": i[4],
            "Stimulus": Stimulus,
            "Perception": perception,
            "Objective": (i[5]),
            "Question": answer}, ignore_index = True)


#Save the logfile
logfile_name = "Cognition and communication exam/logfile/{}_{}_logfile.csv".format(ID, Navn) 
logfile.to_csv(logfile_name)