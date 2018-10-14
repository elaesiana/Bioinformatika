#%%
from Bio.HMM import MarkovModel, Trainer, Utilities
from Bio.Seq import MutableSeq, Seq
from Bio import Alphabet
import random

class state(Alphabet.Alphabet):
    letters = ['1', '2', '3']
class DNA(Alphabet.Alphabet):
    letters = ['A','C','G','T']

model = MarkovModel.MarkovModelBuilder(state(), DNA())
model.allow_all_transitions()

#set probailitas awal secara random
model.set_random_probabilities()

baum_welch = model.get_markov_model()

VERBOSE = 0
def stop_training(log_likelihood_change, n_iterasi):
    if VERBOSE:
        print("ll change: %f" % log_likelihood_change)
    if log_likelihood_change < 0.01:
        return 1
    elif n_iterasi >= 10:
        return 1
    else:
        return 0


def random_generator(num):
    states = MutableSeq('',state())
    for i in range(num):
        states.append(random.choice('123'))
    

    sequence = MutableSeq('',DNA())
    for i in range(num):
        sequence.append(random.choice('ACTG'))
    
    return states.toseq(),sequence.toseq()

n_seq = 30
states = []
sequence = []
seq = []

for i in range(n_seq):
    num = 10
    state2, sequence2 = random_generator(num)
    states.append(state2)
    sequence.append(sequence2)
    seq.append(Trainer.TrainingSequence(sequence[i],states[i]))
    
trainer = Trainer.BaumWelchTrainer(baum_welch)
trained = trainer.train(seq, stop_training)

print('\n\nProbabilitas Transisi: ',trained.transition_prob)
print('\nProbabilitas Emisi: ',trained.emission_prob)

prediction, prob = trained.viterbi(sequence[0], state())
print('\nProbabilitas Prediksi: ', prob)
Utilities.pretty_print_prediction(sequence[0], states[0], prediction)
