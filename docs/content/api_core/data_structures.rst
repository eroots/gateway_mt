data_structures
===============

.. automodule:: gateway_mt.data_structures
   :no-members:

.. _Dataset:

Dataset
-------

.. autoclass:: gateway_mt.data_structures.Dataset

.. _Data:

Data
----

.. autoclass:: gateway_mt.data_structures.Data

.. _Error Floors:

Error Floors
------------

Error floors may be set in a Data object by setting the corresponding value in the error_floors attribute.
Data.error_floors initialized as::

	Data.error_floors = {'Off-Diagonal Impedance': 0.03,
	                     'Diagonal Impedance': 0.03,
	                     'Tipper': 0.03,
	                     'Rho': 0.03,
	                     'Phase': 0.02}

The corresponding entries may be modified as any dictionary and then accesing the apply_error_floor method. For instance, to change the error floor of the ZXY and ZYX components to 7.5%, use::

	data.error_floors['Off-Diagonal Impedance'] = 0.075
	data.apply_error_floor()

.. _RawData:

RawData
-------

.. autoclass:: gateway_mt.data_structures.RawData

.. _Response:

Response
--------

.. autoclass:: gateway_mt.data_structures.Response

.. _Site:

Site
----

.. autoclass:: gateway_mt.data_structures.Site

.. _Model:

Model
-----

.. autoclass:: gateway_mt.data_structures.Model

.. _Phase Tensors:

PhaseTensor
-----------

.. autoclass:: gateway_mt.data_structures.PhaseTensor

.. autoclass:: gateway_mt.data_structures.CART
