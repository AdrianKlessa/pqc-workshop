## Algorithms and helper functions related to post-quantum cryptography

### Disclaimer:

### This code is not safe for real-world usage

*The implementations are for educational purposes*

*They are slow, don't use proper sources of entropy and are surely vulnerable to timing and various other attacks*

The implementation of Kyber used here also doesn't utilize the Number Theoretic Transform (NTT) so it's both slower and not compliant with the NIST specification of the algorithm.

The Kyber code also doesn't verify that the message's error vectors were generated properly.

Useful sources used for implementing Kyber:

* https://cryptopedia.dev/posts/kyber/
* https://pq-crystals.org/kyber/data/kyber-specification-round3-20210131.pdf
* https://cryptographycaffe.sandboxaq.com/posts/kyber-01/
