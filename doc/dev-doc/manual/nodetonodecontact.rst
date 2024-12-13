.. _sect-ntn:

Node-to-node contact
--------------------

Rate-and-state friction laws
````````````````````````````

Dieterich law
'''''''''''''

A "+1" term is added inside the log of velocity to avoid NaN values when the velocity tends toward 0.

.. math::
  \mu = \mu_0 + a \ln \left(1 + \frac{|v|}{v^*}\right) + b\ln \left(\frac{\phi}{\phi^*}\right)


Revised Dieterich law
'''''''''''''''''''''

A "+1" term is added inside the log of velocity to avoid NaN values when the velocity tends toward 0. Note the :math:`b \mu_0` to have a coherent shape with the N-shaped law.

.. math::
  \mu = \mu_0 + a \ln \left(1 + \frac{|v|}{v^*}\right) + b\mu_0\ln \left(\frac{\phi}{\phi^*}\right)

Regularised Dieterich law
'''''''''''''''''''''''''

A "+1" term is added inside the log of velocity to avoid NaN values when the velocity tends toward 0.

.. math::
  \mu = \mu_0 + a \ln \left(1 + \frac{|v|}{v^*}\right) + b\mu_0\ln \left(1 + \frac{\phi}{\phi^*}\right)

N-shaped law
''''''''''''

The N-shaped law features rate-weakening at intermediate velocities and a strengthening behavior at really low and high sliding velocities. The "+1" term in the logarithm of the state leads to the existence of a velocity-strengthening branch at large velocity. 

.. math::
  \mu = \left[1 + b \ln \left(1 + \frac{\phi}{\phi*} \right)\right] \left[\frac{\mu_0}{\sqrt{1 + (v^*/v)^2}} + a \ln \left(\frac{|v|}{v^*}\right)\right]

Revised N-shaped law
''''''''''''''''''''

The N-shaped law features rate-weakening at intermediate velocities and a strengthening behavior at really low and high sliding velocities.
The "+1" term in the logarithm of the velocity and the function :math:`\sqrt{1 + (v^*/v)^2}` ensures that the friction coefficient vanishes for :math:`v=0` and that the behavior at very small velocity is rate-strengthening (this regime is related to the creep of microcontacts). The "+1" term in the logarithm of the state leads to the existence of a velocity-strengthening branch at large velocity. 

.. math::
  \mu = \left[1 + b \ln \left(1 + \frac{\phi}{\phi*} \right)\right] \left[\frac{\mu_0}{\sqrt{1 + (v^*/v)^2}} + a \ln \left(1 + \frac{|v|}{v^*}\right)\right]

