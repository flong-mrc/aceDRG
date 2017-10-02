import os, sys
sys.argv[0] = os.path.join(os.environ['CCP4'], 'bin', 'acedrg')
import acedrg
acedrg.main()
