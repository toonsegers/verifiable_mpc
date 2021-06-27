# Verifiable MPC
The 'Verifiable MPC' Python package implements the verifiable secure multi-party computation (MPC) scheme.

## Verifiable MPC scheme
Electronic voting protocols are prime examples of secure multi-party computation (MPC) that separate input and compute parties. In this setting, there are many input parties that outsource the tallying operation to a set of compute parties. This setting introduces new requirements versus classical MPC protocols where parties are considered both input and compute party. A necessary requirement for voting protocols is public verifiability. While voting protocols specialize in the linear operation of tallying votes, the focus of this work is a scheme that defines *publicly verifiable MPC for general arithmetic circuits*.

The compute parties produce a zero-knowledge proof of correctness of the computation that allows anyone, particularly someone external to the secure computation, to check the correctness of the output, while preserving the privacy properties of the MPC protocol.

Our scheme addresses the challenge of general arithmetic circuits using recent results in zero-knowledge proof systems, particularly compressed Sigma-protocols [https://eprint.iacr.org/2020/152] by Attema and Cramer (AC20), and Bulletproofs [https://eprint.iacr.org/2017/1066] by Bünz, Bootle, Boneh, Poelstra, Wuille and Maxwell (BBB+17). 

Our construction is based on AC20, which reconciles Bulletproofs with Sigma-protocol theory. The construction yields proofs of logarithmic size and does not require a trusted setup, i.e., the setup does not require knowledge of a trapdoor. 

For our implementation we use the MPyC framework: [https://github.com/lschoe/].

## Installation

This implementation depends on MPyC (version 0.74 or above), gmpy2 and Secure Groups [https://github.com/toonsegers/sec_groups/].

Install latest version of MPyC:

	git clone https://github.com/lschoe/mpyc
	cd mpyc
	python setup.py install

Install 'gmpy2':

	pip install gmpy2   				# for Linux (first running `apt install libmpc-dev` may be necessary)
	pip install gmpy2-[version etc].whl	# for Windows, see Gohlke's unofficial binaries [https://www.lfd.uci.edu/~gohlke/pythonlibs/]

Install Secure Groups:

	git clone https://github.com/toonsegers/sec_groups/
	cd sec_groups
	python setup.py install

## Demos

The following demos are included:

* `demo_circuit_builder.py` to use standard Python and automatically construct an arithmetic circuit in memory;
TODO

Run the demos as follows:

	cd demos
	python demo_circuit_builder.py

## Testing

Run the following commands:

	python -m unittest discover .

## Acknowledgements

This work has received funding from the European Union's Horizon 2020 research and innovation program under grant agreements No 780477 (PRIViLEDGE).
