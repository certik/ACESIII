
      integer max_eps
      parameter (max_eps = 2000)
      common /epsilon/epsilon(max_eps), epsilonb(max_eps)
      double precision epsilon, epsilonb

      common /fock/focka(max_eps), fockb(max_eps)
      double precision focka, fockb
