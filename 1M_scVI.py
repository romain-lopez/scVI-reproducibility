print "STARTING" 
import tensorflow as tf
from benchmarking import *
import numpy as np
import time

import scVI
from benchmarking import *
from helper import *

data_path = "/home/ubuntu/single-cell-scVI/data/10x1M/data_small.hdf5"
import h5py

f = h5py.File(data_path)    
expression_train = f["data_train"][:1000000]
expression_test = f["data_test"][:]

log_library_size = np.log(np.sum(expression_train, axis=1))
mean, var = np.mean(log_library_size), np.var(log_library_size)

batch_size = 128
learning_rate = 0.0005
epsilon = 0.01

print "building graph"

tf.reset_default_graph()
expression = tf.placeholder(tf.float32, (None, expression_train.shape[1]), name='x')
kl_scalar = tf.placeholder(tf.float32, (), name='kl_scalar')
optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate, epsilon=epsilon)
training_phase = tf.placeholder(tf.bool, (), name='training_phase')

model = scVI.scVIModel(expression=expression, kl_scale=kl_scalar, \
                         optimize_algo=optimizer, phase=training_phase, \
                          library_size_mean=mean, library_size_var=var, \
                          n_layers=3, n_hidden=256)

# Session creation
sess = tf.Session()
sess.run(tf.global_variables_initializer())

print "\n \n TRAINING"
begin = time.time()

result, rho_early, rho_final = train_model_patience(model, (expression_train, expression_test), sess, 120)

end = time.time()
print "time", begin-end
print result

np.save("rho_early", rho_early)
np.save("rho_final", rho_final)

