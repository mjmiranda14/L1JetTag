#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import tensorflow as tf
import tensorflow 
from keras.models import Model
from tensorflow.keras.layers import Conv1D, Dense, Flatten, Input, GlobalAveragePooling1D

x = inputs = Input(shape=(10, 14))
x = Conv1D(
    filters=50,
    kernel_size=1,
    strides=1,
    activation="relu",
)(x)
x = Conv1D(filters=50, kernel_size=1, activation="relu")(x)

# Tensor([10, 50]) -> Tensor([50])
# x = (x[0,:] + x[1,:] + x[2,:] + ...)/10
x = GlobalAveragePooling1D()(x)

# Tensor([10, 50]) -> Tensor([500])
# x = Flatten()(x)

x = Dense(50, activation="relu")(x)
x = Dense(10, activation="relu")(x)
outputs = Dense(1, activation="sigmoid")(x)

model = Model(inputs=inputs, outputs=outputs)
model.summary()
tensorflow.keras.utils.plot_model(model, show_shapes=True)

import hls4ml
config = hls4ml.utils.config_from_keras_model(model, granularity='name')

config['LayerName']['conv1d']['ReuseFactor'] = 10
config['LayerName']['conv1d_relu']['ReuseFactor'] = 10
config['LayerName']['conv1d_1']['ReuseFactor'] = 10
config['LayerName']['conv1d_1_relu']['ReuseFactor'] = 10

config['LayerName']['dense']['ReuseFactor'] = 10
# config['LayerName']['dense']['Strategy'] = 'Resource'

hls_model = hls4ml.converters.convert_from_keras_model(model,
                                                       hls_config=config,
                                                       output_dir='model_1/hls4ml_prj',
                                                       part='xcu250-figd2104-2L-e',
                                                       default_precision="ap_fixed<8,3>")

hls4ml.utils.plot_model(hls_model, show_shapes=True, show_precision=True, to_file=None)

hls_model.compile() # spits out all the HLS C++ code

hls_model.build(csim=False, synth=True, cosim=False, validation=False, export=False) # tries to synthesize into RTL



