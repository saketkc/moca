import numpy as np
import cPickle
from pymongo import MongoClient
from bson.binary import Binary
import base64

def encode_image(image_path):
    """Base64 encode image to store in mongodb

    Parameters
    ----------
    image_path: str
        Path to image file

    Returns
    -------
    encoded_img: base64encode
        base64 encoded image
    """
    with open(image_path, 'rb') as image_fh:
        encoded_img = base64.b64encode(image_fh.read())
    return encoded_img

def pickle_numpy_array(np_array):
    """Create a pickle instance of numpy.array

    Parameter
    ---------
    np_array: np.array
        numpy array

    Returns
    -------
    pickle_field: pickled instance
        Pickled instance
    """
    return cPickle.dumps(np_array)

def unpickle_numpy_array(pickle_field):
    """Unpickle a pickled instance of numpy.array

    Parameter
    ---------
    pickle_field: pickled instance
        Pickled instance

    Returns
    -------
    np_array: np.array
        numpy array
    """
    return cPickle.loads(pickle_field)

def create_binary_pickle(np_array):
    """Create binary version of pickled array
    used for inserting in mongodb

    Parameter
    ---------
    np_array: np.array
        numpy array

    Returns
    -------
    pickle_field: Binary pickled instance
        Binary Pickled instance
    """
    return Binary(pickle_numpy_array(np_array))
