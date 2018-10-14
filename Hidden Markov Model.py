from Bio.HMM import MarkovModel, Trainer, Utilities
from Bio.Seq import MutableSeq, Seq
from Bio import Alphabet
import random

class state(Alphabet.Alphabet):
    letters = ['1', '2', '3']
class observation(Alphabet.Alphabet):
    letters = ['A','C','G','T']

mm = MarkovModel.MarkovModelBuilder(state(), observation())
mm.allow_all_transitions()

#set probailitas awal secara random
mm.set_random_probabilities()

baum_welch = mm.get_markov_model()

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


num = 200
states = MutableSeq('',state())
for i in range(num):
    states.append(random.choice('123'))
states.toseq()

alphabet = MutableSeq('',observation())
for i in range(num):
    alphabet.append(random.choice('ACTG'))
alphabet.toseq()


seq = Trainer.TrainingSequence(alphabet,states)
trainer = Trainer.BaumWelchTrainer(baum_welch)
trained = trainer.train([seq], stop_training)

print('\n\nProbabilitas Transisi: ',trained.transition_prob)
print('\nProbabilitas Emisi: ',trained.emission_prob)

prediction, prob = trained.viterbi(alphabet, state())
print('\nProbabilitas Prediksi: ', prob)
Utilities.pretty_print_prediction(alphabet, states, prediction)
