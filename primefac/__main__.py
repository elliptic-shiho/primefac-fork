from .primefac import primefac_cli_main
import sys

primefac_cli_main(sys.argv)


"""
Fun examples:
primefac -v 1489576198567193874913874619387459183543154617315437135656
  On my system, the factor race is a bit unpredicatble on this number.
  prb, ecm, p-1, and mpqs all show up reasonably often.
primefac -v 12956921851257164598146425167654345673426523793463
  Z50 = P14 x P17 x P20 =
      24007127617807 x 28050585032291527 x 19240648901716863967.
  p-1 gets the P14 and p+1 gets the rest.
primefac -v 38 ! 1 +
    -->  Z45 = P23 x P23 = 14029308060317546154181 x 37280713718589679646221
  The MPQS (almost always) gets this one.
    Depending on the system running things,
        this can take from 10 seconds to 3 minutes.
"""

"""
main(['p', '-s',
'1489576198567193874913874619387459183543154617315437135656']) only test
"""
