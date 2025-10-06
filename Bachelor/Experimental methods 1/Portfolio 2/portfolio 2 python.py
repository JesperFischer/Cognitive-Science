#loading modules:
from psychopy import visual, event, core, gui
import random
import pandas as pd
# if we use as afterwards we can referer back to pandas as pd.

#pop up window to get information.
pop_up = gui.Dlg(title = "The short story experiment")
pop_up.addField("Participant ID:")
pop_up.addField("Age:")
pop_up.addField("Gender:", choices = ["female", "male"])

#showing the box to the participant
pop_up.show()

#storing the data from the pop_up window in varibles, and making the experiment quit if cancel is clicked.
if pop_up.OK:
    ID = pop_up.data[0]
    age = pop_up.data[1]
    gender = pop_up.data[2]
elif pop_up.Cancel:
    core.quit()

#defining a function for the start and end message (draw, flip eventkeys).
def d_f_ek(x):
    x.draw()
    win.flip()
    event.waitKeys(keyList=["space","q"])

#defining a function that draws the one word from the story and then flips that to the participants.
def drawflip(x):
    x = visual.TextStim(win, text = word)
    x.draw()
    win.flip()

#define a window
win = visual.Window(fullscr = True, color = "black")

#defining the logfile
logfile = pd.DataFrame(columns = ["ID", "Age", "Gender", "Condition", "Word", "Reaction_time"])

#define a stopwatch:
stopwatch = core.Clock()

#the short story look that there is (.split()) in the end so words come one at a time
condfun = "there once was a fox who was friend with an ox they went to the stream to clean their rings and what was to be found but a cat with a rat The cat was playing with the rat so the fox said to the ox perhaps we can play their way So the fox leaped on top of the ox and the ox threw the fox in the way he knew the fox was high in the air When he decided to declare that he was done and did not find this very fun So he hopped off the ox and then they went home ".split()
condcool = "there once was a fox who was friend with an ox they went to the stream to clean their rings and what was to be found but a cat with a rat The cat was playing with the rat so the fox said to the ox perhaps we can play their way So the fox leaped on top of the ox and the ox threw the fox in the way he knew the fox was high in the air When he decided to declare that he was done and did not find this very fun So he hopped off the cow and then they went home".split()

#defining the welcome message and showing it using the function
msg = visual.TextStim(win, text = "Welcome to the experiment. In this experiment you will be presented with a story that you have to read and understand. Use the space bar to get the next word. Press space to continue")
d_f_ek(msg)


#the two conditions defined from the participant ID, if the ID is an even number, the modulus is 0, and if it is an uneven number the modulus is different from 0. The two conditions are either cow or ox, since that is the difference in the stories.
#here the first condition is difined as ox and is "true" if the ID is an even number (modulus == 0)
if int(ID[0]) % 2 == 0:
    condition = "ox"
    
    #looping though the first condition
    for word in condfun:
        
        #using the previously defined function to draw each word to the participant
        drawflip(word)
        
        #resetting the stopwatch
        stopwatch.reset()
        
        #defining the keys one can press.
        key = event.waitKeys(keyList=["space","q"])
        
        #defining reationtime as the time between the word is shown to the participant till the participant presses a key. 
        reaction_time = stopwatch.getTime()
        
        #putting all the data into the logfile
        logfile = logfile.append({
        "ID": ID, 
        "Age": age, 
        "Gender": gender,
        "Condition": condition,
        "Word": word,
        "Reaction_time": reaction_time}, ignore_index = True)
        
        #making it possible to quit the eksperiment by pressing q.
        if key[0] == "q":
            core.quit()

#here the second condition is defined and is if the participant ID is an uneven number (modulus !=(is different from) 0) and is called cow.
elif int(ID[0]) % 2 != 0:
    condition = "cow"
    
    #looping though the text with the cool word and showing each word.
    for word in condcool:
        
        #using the previously defined function to draw each word to the participant
        drawflip(word)
        
        #resetting the stopwatch
        stopwatch.reset()
        
        #defining the keys one can press.
        key = event.waitKeys(keyList=["space","q"])
        
        #defining reationtime as the time between the word is shown to the participant till the participant presses a key. 
        reaction_time = stopwatch.getTime()
        
        #putting all the data into the logfile
        logfile = logfile.append({
        "ID": ID, 
        "Age": age, 
        "Gender": gender,
        "Condition": condition,
        "Word": word,
        "Reaction_time": reaction_time}, ignore_index = True)
        #making it possible to quit the eksperiment by pressing q.
        if key[0] == "q":
            core.quit()


#defining the end message and showing it using the function:
msg2 = visual.TextStim(win, text = "Thank you for your participation. Press space to end")
d_f_ek(msg2)

#putting the data into the logfile. And saving the logfile.
logfile_name = "logfiles_{}.csv".format(ID)
#saving
logfile.to_csv(logfile_name)
