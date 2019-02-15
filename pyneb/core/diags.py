""" @module
This module contains the admitted diagnostics.
Diagnostics are stored in a dictionary. Keys are literal descriptions of the diagnostic 
(e.g. '[NI] 5198/5200') and items are tuples including a label ('N1'), an expression for the line ratio ('L(5198)/L(5200)',
and an expression for the uncertainty on the line ratio as a function of the uncertainty on the individual lines ('RMS([E(5200), E(5198)])')

"""
import numpy as np
import pyneb as pn
from pyneb.utils.misc import int_to_roman, parseAtom, parseAtom2
from pyneb.utils.init import BLEND_LIST

diags_dict = {}

diags_dict['[CIII] 1909/1907'] = ('C3', 'L(1909)/L(1907)', 'RMS([E(1909), E(1907)])')
diags_dict['[NI] 5198/5200'] = ('N1', 'I(3, 1)/I(2, 1)', 'RMS([E(5200), E(5198)])')
diags_dict['[NII] 5755/6548'] = ('N2', 'L(5755)/L(6548)', 'RMS([E(6548), E(5755)])')
diags_dict['[NII] 5755/6584'] = ('N2', 'L(5755)/L(6584)', 'RMS([E(6584), E(5755)])')
diags_dict['[NII] 5755/6584+'] = ('N2', 'L(5755)/(L(6548)+L(6584))', 'RMS([E(6548)*L(6548)/(L(6548)+L(6584)), E(6584)*L(6584)/(L(6584)+L(6548)), E(5755)])')
diags_dict['[NII] 121m/20.5m'] = ('N2', 'L(1214747)/L(2054427)', 'RMS([E(2054427)/E(1214747)])')
diags_dict['[NIII] 1751+/57.4m'] = ('N3', 'B("1751A+")/L(574000)', 'RMS([E(574000), BE("1751A+")])')
#diags_dict[''] = ('N3','L(1749)/L(1752)','RMS([E(1749), E(1752)])')
diags_dict['[OI] 63m/147m'] = ('O1', 'L(632000)/L(1455000)', 'RMS([E(632000), E(1455000)])')
diags_dict['[OI] 5577/6302'] = ('O1', 'L(5577)/L(6300)', 'RMS([E(6300), E(5577)])')
diags_dict['[OI] 5577/6300'] = ('O1', 'L(5577)/L(6300)', 'RMS([E(6300), E(5577)])')
diags_dict['[OI] 5577/6300+'] = ('O1', 'L(5577)/(L(6300)+L(6364))', 'RMS([E(6300)*L(6300)/(L(6300)+L(6364)), E(6364)*L(6364)/(L(6300)+L(6364)), E(5577)])')
diags_dict['[OII] 3726/3729'] = ('O2', 'L(3726)/L(3729)', 'RMS([E(3729), E(3726)])')
diags_dict['[OII] 3727+/7325+c'] =  ('O2', '(B("3727A+"))/(B("7319A+")+B("7330A+"))',
            'RMS([BE("7319A+")*B("7319A+")/(B("7319A+")+B("7330A+")), BE("7330A+")*B("7330A+")/(B("7319A+")+B("7330A+")), BE("3727A+")])')
diags_dict['[OII] 3727+/7325+'] = ('O2', '(L(3726)+L(3729))/(B("7319A+")+B("7330A+"))',
              'RMS([E(3726)*L(3726)/(L(3726)+L(3729)), E(3729)*L(3729)/(L(3726)+L(3729)),BE("7319A+")*B("7319A+")/(B("7319A+")+B("7330A+")),BE("7330A+")*B("7330A+")/(B("7319A+")+B("7330A+"))])')
#diags_dict['[OII] 3727+/7325+'] = ('O2', '(L(3726)+L(3729))/(B("7325A+"))', 'RMS([E(3726)*L(3726)/(L(3726)+L(3729)), E(3729)*L(3729)/(L(3726)+L(3729)),BE("7325A+")])')
diags_dict['[OII] 3727+/7325+b'] = ('O2', '(L(3726)+L(3729))/(I(4,2)+I(4,3)+I(5,2)+I(5,3))',
              'RMS([E(3726)*L(3726)/(L(3726)+L(3729)), E(3729)*L(3729)/(L(3726)+L(3729)),BE("7319A+")*B("7319A+")/(B("7319A+")+B("7330A+")),BE("7330A+")*B("7330A+")/(B("7319A+")+B("7330A+"))])')
diags_dict['OII 4649.13/4089.29'] = ('O2r', "S('4649.13')/S('4089.29')", "RMS([SE('4649.13'), SE('4089.29')])")
diags_dict['OII 4649.13/4661.63'] = ('O2r', "S('4649.13')/S('4661.63')", "RMS([SE('4649.13'), SE('4661.63')])")
diags_dict['[OIII] 4363/5007'] = ('O3', 'L(4363)/L(5007)', 'RMS([E(5007), E(4363)])')
diags_dict['[OIII] 4363/5007+'] = ('O3', 'L(4363)/(L(5007)+L(4959))', 'RMS([E(5007)*L(5007)/(L(5007)+L(4959)), E(4959)*L(4959)/(L(5007)+L(4959)), E(4363)])')
diags_dict['[OIII] 5007/88m'] = ('O3', 'L(5007)/L(883000)', 'RMS([E(883000), E(5007)])')
diags_dict['[OIII] 51m/88m'] = ('O3', 'L(518000)/L(883000)', 'RMS([E(883000), E(518000)])')
diags_dict['[OIII] 1666/5007+'] = ('O3', 'L(1666)/(L(5007)+L(4959))', 'RMS([E(5007)*L(5007)/(L(5007)+L(4959)), E(4959)*L(4959)/(L(5007)+L(4959)), E(1666)])')
diags_dict['[OIII] 1666/5007'] = ('O3', 'L(1666)/L(5007)', 'RMS([E(5007), E(1666)])')
diags_dict['[OIII] 1664+/5007'] = ('O3', '(B("1664A+"))/L(5007)', 'RMS([BE("1664A+"), E(5007)])')
diags_dict['[OIII] 1666/4363'] = ('O3', 'L(1666)/L(4363)', 'RMS([E(4363), E(1666)])')
diags_dict['[OIV] 1401/1405'] = ('O4', 'L(1401)/L(1405)', 'RMS([E(1401), E(1405)])')
diags_dict['[OIV] 1400+/25.9m'] = ('O4', 'B("1400A+")/L(259000)', 'RMS([BE("1400A+"), E(259000)])')
diags_dict['[SII] 6731/6716'] = ('S2', 'L(6731)/L(6716)', 'RMS([E(6716), E(6731)])')
diags_dict['[SII] 4069/4076'] = ('S2', 'L(4069)/L(4076)', 'RMS([E(4069), E(4076)])')
diags_dict['[SII] 4072+/6720+'] = ('S2', '(L(4069)+L(4076))/(L(6716)+L(6731))',
              'RMS([E(6716)*L(6716)/(L(6716)+L(6731)), E(6731)*L(6731)/(L(6716)+L(6731)), E(4069)*L(4069)/(L(4069)+L(4076)), E(4076)*L(4076)/(L(4069)+L(4076))])')
diags_dict['[SIII] 18.7m/33.5m'] = ('S3', 'L(187000)/L(335000)', 'RMS([E(335000), E(187000)])')
diags_dict['[SIII] 6312/18.7m'] = ('S3', 'L(6312)/L(187000)', 'RMS([E(187000), E(6312)])')
diags_dict['[SIII] 9069/18.7m'] = ('S3', 'L(9069)/L(187000)', 'RMS([E(187000), E(9069)])')
diags_dict['[SIII] 6312/9200+'] = ('S3', 'L(6312)/(L(9069)+L(9531))', 'RMS([E(9069)*L(9069)/(L(9069)+L(9531)), E(9531)*L(9531)/(L(9069)+L(9531)), E(6312)])')
diags_dict['[SIII] 6312/9069'] = ('S3', 'L(6312)/L(9069)', 'RMS([E(9069), E(6312)])')
diags_dict['[NeIII] 15.6m/36.0m'] = ('Ne3', 'L(156000)/L(360000)', 'RMS([E(156000), E(360000)])')
diags_dict['[NeIII] 3869/15.6m'] = ('Ne3', 'L(3869)/L(156000)', 'RMS([E(156000), E(3869)])')
diags_dict['[NeIII] 3930+/15.6m'] = ('Ne3', '(L(3869)+L(3968))/L(156000)', 'RMS([E(156000), E(3869)*L(3869)/(L(3869)+L(3968)), E(3968)*L(3968)/(L(3869)+L(3968))])')
diags_dict['[NeIII] 3344/3930+'] = ('Ne3', 'L(3343)/(L(3869)+L(3968))', 'RMS([E(3869)*L(3869)/(L(3869)+L(3968)), E(3968)*L(3968)/(L(3869)+L(3968)), E(3343)])')
diags_dict['[NeIII] 3343/3930+'] = ('Ne3', 'L(3343)/(L(3869)+L(3968))', 'RMS([E(3869)*L(3869)/(L(3869)+L(3968)), E(3968)*L(3968)/(L(3869)+L(3968)), E(3343)])')
diags_dict['[NeV] 2973/3370+'] = ('Ne5', 'L(2973)/(L(3426)+L(3346))', 'RMS([E(3426)*L(3426)/(L(3426)+L(3346)), E(3346)*L(3346)/(L(3426)+L(3346)), E(2973)])')
diags_dict['[NeV] 1575/3426'] = ('Ne5', 'L(1575)/L(3426)', 'RMS([E(1575), E(3426)])')
diags_dict['[NeV] 14.3m/24.2m'] = ('Ne5', 'L(143000)/L(242000)', 'RMS([E(143000), E(242000)])')
diags_dict['[NeV] 3426/24.2m'] = ('Ne5', 'L(3426)/L(242000)', 'RMS([E(3426), E(242000)])')
diags_dict['[ClIII] 5538/5518'] = ('Cl3', 'L(5538)/L(5518)', 'RMS([E(5518), E(5538)])')
diags_dict['[ClIV] 5323/7531'] = ('Cl4', 'L(5323)/L(7531)', 'RMS([E(7531), E(5323)])')
diags_dict['[ClIV] 5323/7700+'] = ('Cl4', 'L(5323)/(L(7531)+L(8046))', 'RMS([E(7531)*L(7531)/(L(7531)+L(8046)), E(8046)*L(8046)/(L(7531)+L(8046)), E(5323)])')
diags_dict['[ArIII] (7751+7136)/9m'] = ('Ar3', '(L(7751)+L(7136))/L(90000)', 'RMS([E(90000), E(7751)*L(7751)/(L(7751)+L(7136)), E(7136)*L(7136)/(L(7751)+L(7136))])')
diags_dict['[ArIII] 7136/9m'] = ('Ar3', 'L(7136)/L(90000)', 'RMS([E(90000), E(7136)])')
diags_dict['[ArIII] 5192/7300+'] = ('Ar3', 'L(5192)/(L(7751)+L(7136))', 'RMS([E(7751)*L(7751)/(L(7751)+L(7136)), E(7136)*L(7136)/(L(7751)+L(7136)), E(5192)])')
diags_dict['[ArIII] 5192/7136'] = ('Ar3', 'L(5192)/L(7136)', 'RMS([E(7136), E(5192)])')
diags_dict['[ArIII] 9.0m/21.8m'] = ('Ar3', 'L(89897)/L(218000)', 'RMS([E(89897), E(218000)])')
diags_dict['[ArIV] 4740/4711'] = ('Ar4', 'L(4740)/L(4711)', 'RMS([E(4711), E(4740)])')
diags_dict['[ArIV] 2860+/4720+'] = ('Ar4', '(L(2854)+L(2868))/(L(4711)+L(4740))',
              'RMS([E(4711)*L(4711)/(L(4711)+L(4740)), E(4740)*L(4740)/(L(4711)+L(4740)), E(2854)*L(2854)/(L(2854)+L(2868)), E(2868)*L(2854)/(L(2854)+L(2868))])')
diags_dict['[ArIV] 7230+/4720+'] = ('Ar4', '(L(7170)+L(7263))/(L(4711)+L(4740))',
              'RMS([E(4711)*L(4711)/(L(4711)+L(4740)), E(4740)*L(4740)/(L(4711)+L(4740)), E(7170)*L(7170)/(L(7170)+L(7263)), E(7263)*L(7263)/(L(7170)+L(7263))])')
diags_dict['[ArV] 4626/6600+'] = ('Ar5', 'L(4626)/(L(6435)+L(7005))', 'RMS([E(6435)*L(6435)/(L(6435)+L(7005)), E(7005)*L(7005)/(L(6435)+L(7005)), E(4626)])')

diags_dict['[FeIII] 5270/4987'] = ('Fe3', 'L(5270)/L(4987)', 'RMS([E(5270), E(4987)])')
diags_dict['[FeIII] 5270/4925'] = ('Fe3', 'L(5270)/L(4925)', 'RMS([E(5270), E(4925)])')
diags_dict['[FeIII] 5270/4881'] = ('Fe3', 'L(5270)/L(4881)', 'RMS([E(5270), E(4881)])')
diags_dict['[FeIII] 5270/5011'] = ('Fe3', 'L(5270)/L(5011)', 'RMS([E(5270), E(5011)])')
diags_dict['[FeIII] 5270/4931'] = ('Fe3', 'L(5270)/L(4931)', 'RMS([E(5270), E(4931)])')
diags_dict['[FeIII] 5270/4659'] = ('Fe3', 'L(5270)/L(4659)', 'RMS([E(5270), E(4659)])')
diags_dict['[FeIII] 5270/4701'] = ('Fe3', 'L(5270)/L(4701)', 'RMS([E(5270), E(4701)])')
diags_dict['[FeIII] 5270/4734'] = ('Fe3', 'L(5270)/L(4734)', 'RMS([E(5270), E(4734)])')
diags_dict['[FeIII] 5270/4009'] = ('Fe3', 'L(5270)/L(4009)', 'RMS([E(5270), E(4009)])')
diags_dict['[FeIII] 4987/4925'] = ('Fe3', 'L(4987)/L(4925)', 'RMS([E(4987), E(4925)])')
diags_dict['[FeIII] 4987/4881'] = ('Fe3', 'L(4987)/L(4881)', 'RMS([E(4987), E(4881)])')
diags_dict['[FeIII] 4987/5011'] = ('Fe3', 'L(4987)/L(5011)', 'RMS([E(4987), E(5011)])')
diags_dict['[FeIII] 4987/4931'] = ('Fe3', 'L(4987)/L(4931)', 'RMS([E(4987), E(4931)])')
diags_dict['[FeIII] 4987/4659'] = ('Fe3', 'L(4987)/L(4659)', 'RMS([E(4987), E(4659)])')
diags_dict['[FeIII] 4987/4701'] = ('Fe3', 'L(4987)/L(4701)', 'RMS([E(4987), E(4701)])')
diags_dict['[FeIII] 4987/4734'] = ('Fe3', 'L(4987)/L(4734)', 'RMS([E(4987), E(4734)])')
diags_dict['[FeIII] 4987/4009'] = ('Fe3', 'L(4987)/L(4009)', 'RMS([E(4987), E(4009)])')
diags_dict['[FeIII] 4925/4881'] = ('Fe3', 'L(4925)/L(4881)', 'RMS([E(4925), E(4881)])')
diags_dict['[FeIII] 4925/5011'] = ('Fe3', 'L(4925)/L(5011)', 'RMS([E(4925), E(5011)])')
diags_dict['[FeIII] 4925/4931'] = ('Fe3', 'L(4925)/L(4931)', 'RMS([E(4925), E(4931)])')
diags_dict['[FeIII] 4925/4659'] = ('Fe3', 'L(4925)/L(4659)', 'RMS([E(4925), E(4659)])')
diags_dict['[FeIII] 4925/4701'] = ('Fe3', 'L(4925)/L(4701)', 'RMS([E(4925), E(4701)])')
diags_dict['[FeIII] 4925/4734'] = ('Fe3', 'L(4925)/L(4734)', 'RMS([E(4925), E(4734)])')
diags_dict['[FeIII] 4925/4009'] = ('Fe3', 'L(4925)/L(4009)', 'RMS([E(4925), E(4009)])')
diags_dict['[FeIII] 4881/5011'] = ('Fe3', 'L(4881)/L(5011)', 'RMS([E(4881), E(5011)])')
diags_dict['[FeIII] 4881/4931'] = ('Fe3', 'L(4881)/L(4931)', 'RMS([E(4881), E(4931)])')
diags_dict['[FeIII] 4881/4659'] = ('Fe3', 'L(4881)/L(4659)', 'RMS([E(4881), E(4659)])')
diags_dict['[FeIII] 4881/4701'] = ('Fe3', 'L(4881)/L(4701)', 'RMS([E(4881), E(4701)])')
diags_dict['[FeIII] 4881/4734'] = ('Fe3', 'L(4881)/L(4734)', 'RMS([E(4881), E(4734)])')
diags_dict['[FeIII] 4881/4009'] = ('Fe3', 'L(4881)/L(4009)', 'RMS([E(4881), E(4009)])')
diags_dict['[FeIII] 5011/4931'] = ('Fe3', 'L(5011)/L(4931)', 'RMS([E(5011), E(4931)])')
diags_dict['[FeIII] 5011/4659'] = ('Fe3', 'L(5011)/L(4659)', 'RMS([E(5011), E(4659)])')
diags_dict['[FeIII] 5011/4701'] = ('Fe3', 'L(5011)/L(4701)', 'RMS([E(5011), E(4701)])')
diags_dict['[FeIII] 5011/4734'] = ('Fe3', 'L(5011)/L(4734)', 'RMS([E(5011), E(4734)])')
diags_dict['[FeIII] 5011/4009'] = ('Fe3', 'L(5011)/L(4009)', 'RMS([E(5011), E(4009)])')
diags_dict['[FeIII] 4931/4659'] = ('Fe3', 'L(4931)/L(4659)', 'RMS([E(4931), E(4659)])')
diags_dict['[FeIII] 4931/4701'] = ('Fe3', 'L(4931)/L(4701)', 'RMS([E(4931), E(4701)])')
diags_dict['[FeIII] 4931/4734'] = ('Fe3', 'L(4931)/L(4734)', 'RMS([E(4931), E(4734)])')
diags_dict['[FeIII] 4931/4009'] = ('Fe3', 'L(4931)/L(4009)', 'RMS([E(4931), E(4009)])')
diags_dict['[FeIII] 4659/4701'] = ('Fe3', 'L(4659)/L(4701)', 'RMS([E(4659), E(4701)])')
diags_dict['[FeIII] 4659/4734'] = ('Fe3', 'L(4659)/L(4734)', 'RMS([E(4659), E(4734)])')
diags_dict['[FeIII] 4659/4009'] = ('Fe3', 'L(4659)/L(4009)', 'RMS([E(4659), E(4009)])')
diags_dict['[FeIII] 4701/4734'] = ('Fe3', 'L(4701)/L(4734)', 'RMS([E(4701), E(4734)])')
diags_dict['[FeIII] 4701/4009'] = ('Fe3', 'L(4701)/L(4009)', 'RMS([E(4701), E(4009)])')
diags_dict['[FeIII] 4734/4009'] = ('Fe3', 'L(4734)/L(4009)', 'RMS([E(4734), E(4009)])')
# To be checked:
# 4989 is to be changed into 4987, but what about 4987?
# which is the relationship between the 3 lines in the numerator?
# is the error right?
# diags_dict['[FeIII] 4659+/4987+'] = ('Fe3', '(L(4659)+L(4734)+L(4009))/(L(4987)+L(4989))', 'RMS([E(4659), E(4987)])')
diags_dict['[NiIII] 6000/7890'] = ('Ni3', 'L(6000)/L(7890)', 'RMS([E(6000), E(7890)])')
diags_dict['[NiIII] 6401/7890'] = ('Ni3', 'L(6401)/L(7890)', 'RMS([E(6401), E(7890)])')
diags_dict['[NiIII] 6534/7890'] = ('Ni3', 'L(6534)/L(7890)', 'RMS([E(6534), E(7890)])')
diags_dict['[NiIII] 6682/7890'] = ('Ni3', 'L(6682)/L(7890)', 'RMS([E(6682), E(7890)])')
diags_dict['[NiIII] 6797/7890'] = ('Ni3', 'L(6797)/L(7890)', 'RMS([E(6797), E(7890)])')
diags_dict['[NiIII] 6946/7890'] = ('Ni3', 'L(6946)/L(7890)', 'RMS([E(6946), E(7890)])')
diags_dict['[NiIII] 6000/8500'] = ('Ni3', 'L(6000)/L(8500)', 'RMS([E(6000), E(8500)])')
diags_dict['[NiIII] 6401/8500'] = ('Ni3', 'L(6401)/L(8500)', 'RMS([E(6401), E(8500)])')
diags_dict['[NiIII] 6534/8500'] = ('Ni3', 'L(6534)/L(8500)', 'RMS([E(6534), E(8500)])')
diags_dict['[NiIII] 6682/8500'] = ('Ni3', 'L(6682)/L(8500)', 'RMS([E(6682), E(8500)])')
diags_dict['[NiIII] 6797/8500'] = ('Ni3', 'L(6797)/L(8500)', 'RMS([E(6797), E(8500)])')
diags_dict['[NiIII] 6946/8500'] = ('Ni3', 'L(6946)/L(8500)', 'RMS([E(6946), E(8500)])')
diags_dict['[NiIII] 6000/7125'] = ('Ni3', 'L(6000)/L(7125)', 'RMS([E(6000), E(7125)])')
diags_dict['[NiIII] 6401/7125'] = ('Ni3', 'L(6401)/L(7125)', 'RMS([E(6401), E(7125)])')
diags_dict['[NiIII] 6534/7125'] = ('Ni3', 'L(6534)/L(7125)', 'RMS([E(6534), E(7125)])')
diags_dict['[NiIII] 6682/7125'] = ('Ni3', 'L(6682)/L(7125)', 'RMS([E(6682), E(7125)])')
diags_dict['[NiIII] 6797/7125'] = ('Ni3', 'L(6797)/L(7125)', 'RMS([E(6797), E(7125)])')
diags_dict['[NiIII] 6946/7125'] = ('Ni3', 'L(6946)/L(7125)', 'RMS([E(6946), E(7125)])')
diags_dict['[NiIII] 6000/6401'] = ('Ni3', 'L(6000)/L(6401)', 'RMS([E(6000), E(6401)])')
diags_dict['[NiIII] 6534/6401'] = ('Ni3', 'L(6534)/L(6401)', 'RMS([E(6534), E(6401)])')
diags_dict['[NiIII] 6946/6401'] = ('Ni3', 'L(6946)/L(6401)', 'RMS([E(6946), E(6401)])')
diags_dict['[NiIII] 6000/6797'] = ('Ni3', 'L(6000)/L(6797)', 'RMS([E(6000), E(6797)])')
diags_dict['[NiIII] 6534/6797'] = ('Ni3', 'L(6534)/L(6797)', 'RMS([E(6534), E(6797)])')
diags_dict['[NiIII] 6946/6797'] = ('Ni3', 'L(6946)/L(6797)', 'RMS([E(6946), E(6797)])')
diags_dict['[NiIII] 6000/6682'] = ('Ni3', 'L(6000)/L(6682)', 'RMS([E(6000), E(6682)])')
diags_dict['[NiIII] 6534/6682'] = ('Ni3', 'L(6534)/L(6682)', 'RMS([E(6534), E(6682)])')
diags_dict['[NiIII] 6946/6682'] = ('Ni3', 'L(6946)/L(6682)', 'RMS([E(6946), E(6682)])')

class Diagnostics(object):
    """
    Diagnostics is the class used to manage the diagnostics and to computed physical conditions 
        (electron temperatures and densities) from them. 
    It is also the class that plots the diagnostic Te-Ne diagrams.

    """    
    def __init__(self, addAll=False, OmegaInterp='Cheb', NLevels=None):
        """
        Diagnostics constructor
        
        Parameters:
            - addAll        switch to include all defined diagnostics (default = False)
            - OmegaInterp   parameter sent to Atom, default is 'Cheb', other can be 'linear'
            
        """
        self.log_ = pn.log_ 
        self.calling = 'Diagnostics'
        ##            
        # @var diags
        # The dictionary containing the diagnostics
        self.diags = {}
        ##
        # @var atomDict
        # The dictionary containing the atoms used for the diagnostics        
        self.atomDict = {}
        self.OmegaInterp = OmegaInterp
        self.NLevels = NLevels
        if addAll:
            self.addAll()


    def getDiagFromLabel(self, label):
        """
        Return the definition of a diagnostic (the 3 or 4 elements tuple)
        
        Usage:
            diags.getDiagFromLabel('[NII] 5755/6548')
        
        Parameters:
            -label a diagnostic label 
            
        """
        if label in self.diags:
            return self.diags[label]
        else:
            self.log_.warn('{0} not in diagnostic'.format(label), calling=self.calling)
            return None

        
    def getDiags(self):
        """
        Return the definitions (tuples) of all the diagnostics defined in the Object.
        No parameter.
        
        """
        return [self.getDiagFromLabel(label) for label in self.diags]

            
    def getDiagLabels(self):
        """
        Return the labels of all the diagnostics defined in the Object
        No parameter.
        
        """
        return self.diags.keys()    

    
    def getAllDiags(self):
        """
        Return the definitions (tuples) of all the possible diagnostics.
        No parameter.
        
        """
        return diags_dict

            
    def getAllDiagLabels(self):
        """
        Return the labels of all the possible diagnostics.
        No parameter.
        
        """
        return diags_dict.keys()    

    
    def getUniqueAtoms(self):
        """
        Return a numpy.ndarray of the ions needed by the diagnostics. Unique is applied to the list before returning.
        No parameter.
        
        """
        return np.unique([self.diags[d][0] for d in self.diags if self.diags[d] is not None])
    
    
    def addDiag(self, label=None, diag_tuple=None, atom=None, wave_range=None):
        """
        Add diagnostics to the list of available diagnostics. It can either be one of the built-in diagnostics,
        a new, user-defined one, or a subset of the built-in diagnostics corresponding to a given atom or wavelength.
        
        Parameters:
            - label        a string or a list of strings describing the diagnostic, e.g. '[OIII] 4363/5007'. 
                           If it is not a key of diags_dict (a diagnostic define by PyNeb), diag_tuple must also be specified
            - diag_tuple   a 3 elements tuple containing:
                           + the atom, e.g. 'Ar5'
                           + the algebraic description of the diagnostic, in terms of line wavelengths or blends or levels, 
                             e.g. '(L(6435)+L(7006))/L(4626)'
                           + the algebraic description of the error, e.g. 
                             'RMS([E(6435)*L(6435)/(L(6435)+L(7006)), E(7006)*L(7006)/(L(6435)+L(7006)), E(4626)])'
            - atom         the selected atom, e.g. 'O3'
            - wave_range   the selected wavelength range
            
        Usage:
        diags.addDiag('[OIII] 4363/5007')
        diags.addDiag('[OIII] 5007/51m', ('O3', 'L(5007)/L(51800)', 'RMS([E(51800), E(5007)])'))
        diags.addDiag(atom='O3')
        diags.addDiag(wave_range=[4000, 6000])

        """
        if type(label) is list:
            for lab in label:
                self.addDiag(lab)
        elif label in self.diags:
            self.log_.warn('{0} already in diagnostic'.format(label), calling=self.calling)
        elif label in diags_dict:
            self.log_.message('Adding diag {0}'.format(label), calling=self.calling)
            self.diags[label] = diags_dict[label]
            atom = diags_dict[label][0]
        elif type(diag_tuple) is tuple:
            if len(diag_tuple) == 3:    
                self.diags[label] = diag_tuple
                atom = diag_tuple[0]
            else:
                self.log_.error('{0} is not in the list of diagnostics. The parameter diag_tuple must be a 3-elements tuple describing the diagnostic'.format(label), calling=self.calling + '.addDiag')
        elif atom is not None:
            for label in diags_dict:
                if diags_dict[label][0] == atom:
                    self.addDiag(label)
#                    if atom not in self.atomDict:
#                        self.atomDict[atom] = pn.Atom(parseAtom(atom)[0], parseAtom(atom)[1])
        elif wave_range is not None:
# To be done: add the possibility to use several independent ranges. The following bit does the trick, but
# only if all the wavelengths of a diagnostic are included in the same range
#            if type(wave_range[0]) is list:
#                for subrange in wave_range:
#                    self.addDiag(label, diag_tuple, atom, subrange)
                 
            for label in diags_dict:
                expr = diags_dict[label][1]
                waves = []
                wave = ''
                in_wave = False
                for i in expr:
                    if i.isdigit():
                        wave = wave + i
                        in_wave = True
                    elif (i == 'm'):
                        if in_wave == True:
                            waves.append(int(wave) * 10000.)
                            in_wave = False
                            wave = ''
                    else:
                        if in_wave == True:
                            waves.append(int(wave))
                            in_wave = False
                            wave = ''
                if (min(waves) > min(wave_range)) and (max(waves) < max(wave_range)):
                    self.addDiag(label)
        else:
            self.log_.error('Bad syntax. You must either give the label of an existing diagnostic, or the label and the tuple of a new one,' + 
                            'or an ion, or a wave range. label={0}, diag_tuple={1}, atom={2}, wave_range={3}'.format(label, diag_tuple, atom, wave_range),
                            calling=self.calling + '.addDiag')
        if atom not in self.atomDict and (type(label) is not list):
            self.atomDict[atom] = pn.Atom(parseAtom(atom)[0], parseAtom(atom)[1], NLevels=self.NLevels)
               
               
    def addAll(self):
        """
        Insert all the possible diagnostics in the Object (no parameters)
        
        """
        for label in diags_dict:
            self.addDiag(label)
            
            
    def delDiag(self, label=None):
        """
        Remove a diagnostic, based on its label.
        
        Parameter:
            - label  the label of the diagnostic ('all': removes all the diagnostics)

        """
        if label == 'all':
            for item in self.diags:
                self.delDiag(item)
        elif label in self.diags:
            del self.diags[label]
        else:
            self.log_.warn('{0} not in diagnostic, cannot be removed'.format(label), calling=self.calling)
    
    
    def addDiagsFromObs(self, obs):
        """
        Add all the possible diagnostics that can be computed from an Observation object
        
        Usage:
            diags.addDiagsFromObs(obs)
            
        Parameter:
            obs     an Observation object
        """
        
        if not isinstance(obs, pn.Observation):
            pn.log_.error('The argument must be an Observation object', calling=self.calling + 'addDiagsFromObs')
        old_level = pn.log_.level
        def I(i, j):
            wave = atom.wave_Ang[i - 1, j - 1]
            corrIntens = obs.getLine(sym, spec, wave).corrIntens
            return corrIntens
        def L(wave):
            corrIntens = obs.getLine(sym, spec, wave).corrIntens
            return corrIntens
        def B(label):
            full_label = atom + '_' + label
            corrIntens = obs.getLine(label=full_label).corrIntens
            return corrIntens
        def S(label):
            full_label = atom + '_' + label + 'A'
            corrIntens = obs.getLine(label=full_label).corrIntens
            return corrIntens
            
        for label in diags_dict:
            atom, diag_expression, error = diags_dict[label]
            sym, spec, rec = parseAtom2(atom)
            if label == '[OII] 3727+/7325+c':
                try:
                    diag_value = eval(diag_expression)
                except Exception as ex:
                    pn.log_.level = old_level
                    pn.log_.debug('Diag not valid {} {}'.format(label, diag_expression))
            try:
                pn.log_.level = 1
                diag_value = eval(diag_expression)
                pn.log_.level = old_level
                if atom not in self.atomDict:
                    if rec == 'r':
                        self.atomDict[atom] = pn.RecAtom(atom=sym+spec)
                    else:
                        self.atomDict[atom] = pn.Atom(atom=atom, NLevels=self.NLevels)
                self.addDiag(label)
            except Exception as ex:
                pn.log_.level = old_level
                pn.log_.debug('Diag not valid {} {}'.format(label, diag_expression))
                
    
    def setAtoms(self, atom_dic):
        """
        Define the dictionary containing the atoms used for the diagnostics.
        A dictionary of atom instantiations refereed by atom strings, for example:
        {'O3': pn.Atom('O', 3)}
        
        """
        if type(atom_dic) != type({}):
            self.log_.error('the parameter must be a dictionary.', calling=self.calling + '.setAtoms')
            return None
        for atom in atom_dic:
            if not isinstance(atom_dic[atom], pn.Atom) and not isinstance(atom_dic[atom], pn.RecAtom):
                self.log_.error('the parameter must be a dictionary of Atom.', calling=self.calling + '.setAtoms')
                return None
        self.atomDict = atom_dic
    
    
    def addClabel(self, label, clabel):
        """
        Add an alternative label to a diagnostic that can be used when plotting diagnostic diagrams.
        
        """
        if label in self.diags:
            self.diags[label] = (self.diags[label][0], self.diags[label][1], self.diags[label][2], clabel)
        else:
            pn.log_.warn('Try to add clabel in undefined label {0}'.format(label), calling=self.calling)
    
    
    def plot(self, emis_grids, obs, quad=True, i_obs=None, alpha=0.3, ax=None, error_band=True):
        """
        PLotting tool to generate Te-Ne diagrams.
        
        Usage:
            diags.plot(emisgrids, obs, i_obs=3)
            
        Parameters:
            - emis_grids    A dictionary of EmisGrid objects refereed by their atom strings (e.g. 'O3')
                            This can for example be the output of pn.getEmisGridDict()
            - obs           A pn.Observation object that is supposed to contain the line intensities
                            used for the plot (corrected intensities).
            - quad          If True (default) the sum of the error is quadratic,otherwise is linear.
            - i_obs         reference for the observation to be plotted, in case there is more than one
                            in the obs object
            - alpha         Transparency for the error bands in the plot
            - error_band    Boolean: plot [default] an error band
            
        """
        if not pn.config.INSTALLED['plt']: 
            pn.log_.error('Matplotlib not available, no plot', calling=self.calling + '.plot')
            return None
        else:
            import matplotlib.pyplot as plt
        if type(emis_grids) != type({}):
            self.log_.error('the first parameter must be a dictionary', calling=self.calling + '.plot')
            return None
        for em in emis_grids:
            if not isinstance(emis_grids[em], pn.EmisGrid):
                self.log_.error('the first parameter must be a dictionary of EmisGrid', calling=self.calling + '.plot')
                return None
        if not isinstance(obs, pn.Observation):
            self.log_.error('the second parameter must be an Observation', calling=self.calling + '.plot')
            return None
        if (i_obs is None) and (obs.n_obs != 1):
            self.log_.error('i_obs must be specified when obs is multiple. try i_obs=0', calling=self.calling)
            return None
        if ax is None:
            f, ax = plt.subplots()
        else:
            f = plt.gcf()
        X = np.log10(emis_grids[list(emis_grids.keys())[0]].den2D)
        Y = emis_grids[list(emis_grids.keys())[0]].tem2D
        for label in self.diags:
            self.log_.message('plotting {0}'.format(label), calling=self.calling)
            diag = self.diags[label]
            atom = diag[0]
            def I(i, j):
                return emis_grids[atom].getGrid(lev_i=i, lev_j=j)
            def L(wave):
                return emis_grids[atom].getGrid(wave=wave)
            def S(label):
                return emis_grids[atom].getGrid(label=label)
            def B(label, I=I, L=L):
                full_label = atom + '_' + label
                if full_label in BLEND_LIST:
                    to_eval = BLEND_LIST[full_label]
                else:
                    self.log_.warn('{0} not in BLEND_LIST'.format(full_label), calling=self.calling)
                    return None
                return eval(to_eval)
            diag_map = eval(diag[1])
            try:
                diag_map = eval(diag[1])
            except:
                diag_map = None
                self.log_.warn('diag {0} {1} not used'.format(diag[0], diag[1]), calling=self.calling)
            if diag_map is not None:
                sym, spec, rec = parseAtom2(atom)
                def I(i, j):
                    wave = emis_grids[atom].atom.wave_Ang[i - 1, j - 1]
                    corrIntens = obs.getLine(sym, spec, wave).corrIntens
                    if i_obs is None:
                        return corrIntens
                    else:
                        return corrIntens[i_obs]
                def L(wave):
                    corrIntens = obs.getLine(sym, spec, wave).corrIntens
                    if i_obs is None:
                        return corrIntens
                    else:
                        return corrIntens[i_obs]
                def S(label):
                    full_label = atom + '_' + label + 'A'
                    corrIntens = obs.getLine(label=full_label).corrIntens
                    if i_obs is None:
                        return corrIntens
                    else:
                        return corrIntens[i_obs]                    
                def B(label):
                    full_label = atom + '_' + label
                    corrIntens = obs.getLine(label=full_label).corrIntens
                    if i_obs is None:
                        return corrIntens
                    else:
                        return corrIntens[i_obs]
                try:
                    ee = eval(diag[1])
                    if type(ee) == np.float64:
                        diag_value = ee
                    else:
                        diag_value = ee[0]
                except:
                    #print(ee, type(ee))
                    diag_value = None
                    self.log_.warn('A line in diagnostic {0} of {1}{2} is missing'.format(diag[1], sym, spec),
                                   calling=self.calling)
                if diag_value is not None and not np.isnan(diag_value) and not np.isinf(diag_value):
                    def E(wave):
                        err = obs.getLine(sym, spec, wave).corrError
                        if i_obs is None:
                            return err
                        else:
                            return err[i_obs]
                    def BE(label):
                        full_label = atom + '_' + label
                        err = obs.getLine(label=full_label).corrError
                        if i_obs is None:
                            return err
                        else:
                            return err[i_obs]
                    def SE(label):
                        full_label = atom + '_' + label + 'A'
                        err = obs.getLine(label=full_label).corrError
                        if i_obs is None:
                            return err
                        else:
                            return err[i_obs]                        
                    if quad is True:
                        RMS = lambda err: np.sqrt((np.asarray(err) ** 2.).sum())
                    else:
                        RMS = lambda err: (np.asarray(err)).sum() 
                    tol_value = eval(diag[2])
                    col_dic = {'C':'cyan', 'N':'blue', 'O':'green', 'Ne':'magenta',
                               'Ar':'red', 'Cl':'magenta', 'S':'black', 'Fe':'blue'}
                    if sym in col_dic:
                        col = col_dic[sym]
                    else:
                        col = 'black'
                    style_dic = {'1':'-', '2':'--', '3':':', '4':'-.', '5':'-', '6':'--'}
                    style = style_dic[spec]
                    if tol_value > 0. and error_band:
                        levels = [(1 - tol_value) * diag_value, (1 + tol_value) * diag_value]
                        if levels[0] < levels[1]:
                            #pn.log_.debug('{} levels {}'.format(label, levels), calling=self.calling)
                            CS = ax.contourf(X, Y, diag_map, levels=levels, alpha=alpha, colors=col)
                    CS = ax.contour(X, Y, diag_map, levels=[diag_value], colors=col, linestyles=style)
                    try:
                        ax.set_xlabel(r'log(n$_{\rm e}$) [cm$^{-3}$]')
                        ax.set_ylabel(r'T$_{\rm e}$ [K]')
                        if len(diag) >= 4:
                            fmt = diag[3]
                        else:
                            fmt = '[{0}{1}]'.format(sym, int_to_roman(int(spec)))
                        ax.clabel(CS, inline=True, fmt=fmt, fontsize=15, colors=col)
                        if type(diag_value) is np.ndarray:
                            diag_value = diag_value[0]
                        self.log_.message('plotted {0}: {1} = {2:.2} with error of {3:.2} %'.format(fmt, label, diag_value, tol_value * 100),
                                          calling=self.calling)
                    except:
                        self.log_.message('NOT plotted {0} {1}'.format(fmt, label),
                                          calling=self.calling)
                        
        
        
    def getCrossTemDen(self, diag_tem, diag_den, value_tem=None, value_den=None, obs=None, i_obs=None,
                       guess_tem=10000, tol_tem=1., tol_den=1., max_iter=5, maxError=1e-3,
                       start_tem= -1, end_tem= -1, start_den= -1, end_den= -1):
        """
        Cross-converge the temperature and density derived from two sensitive line ratios, by inputting the quantity 
        derived with one line ratio into the other and then iterating.
        The temperature- and density-sensitive ratios can be input directly or as an Observation object
    
        Parameters:
    
        - diag_tem   temperature-sensitive diagnostic line ratio
        - diag_den   density-sensitive diagnostic line ratio
        - value_tem  value of the temperature-sensitive diagnostic
        - value_den  value of the density-sensitive diagnostic
        - obs        np.Observation object. Values for observed temperature and density diagnostics are
                        taken from it if value_tem and value_den are None
        - i_obs      index of the observations to be used from obs. 
                        If None, all the observations are considered.
                        May be a boolean array
        - guess_tem  temperature assumed in the first iteration, in K
        - tol_tem    tolerance of the temperature result, in %
        - tol_den    tolerance of the density result, in %
        - max_iter   maximum number of iterations to be performed, after which the function will throw a result
        - maxError   maximum error in the calls to getTemDen, in %
        - start_tem, end_tem  lower and upper limit of the explored temperature range 
        - start_den, end_den  lower and upper limit of the explored density range 
    
        Example:
            tem, den = diags.getCrossTemDen('[OIII] 4363/5007', '[SII] 6731/6716', 0.0050, 1.0, 
                        guess_tem=10000, tol_tem = 1., tol_den = 1., max_iter = 5)

        """
        # TODO:
        # Define a lower/upper density and temperature 
        # to be used in case the diag is pointing to values out of the boundaires
        if diag_tem not in self.diags:
            self.addDiag(diag_tem)
        atom_tem = self.diags[diag_tem][0]
        elem_tem, spec_tem = parseAtom(atom_tem)
        if atom_tem not in self.atomDict:
            self.atomDict[atom_tem] = pn.Atom(elem_tem, spec_tem, self.OmegaInterp, NLevels=self.NLevels)
        atom_tem = self.atomDict[atom_tem]
        if diag_den not in self.diags:
            self.addDiag(diag_den)
        atom_den = self.diags[diag_den][0]
        elem_den, spec_den = parseAtom(self.diags[diag_den][0])
        if (atom_den) not in self.atomDict:
            self.atomDict[atom_den] = pn.Atom(elem_den, spec_den, self.OmegaInterp, NLevels=self.NLevels)
        atom_den = self.atomDict[atom_den]
        eval_tem = self.diags[diag_tem][1]
        eval_den = self.diags[diag_den][1]
        calling = 'Diag.getCrossTemDen %s %s' % (diag_tem, diag_den)
        if value_tem is None:
            def L(wave):
                if i_obs is None:
                    return obs.getLine(elem_tem, spec_tem, wave).corrIntens
                else:
                    return obs.getLine(elem_tem, spec_tem, wave).corrIntens[i_obs]
            def I(i, j):
                wave = atom_tem.wave_Ang[i - 1, j - 1]
                if i_obs is None:
                    return obs.getLine(elem_tem, spec_tem, wave).corrIntens
                else:
                    return obs.getLine(elem_tem, spec_tem, wave).corrIntens[i_obs]
            def B(label, I=I, L=L):
                full_label = elem_tem + spec_tem + '_' + label
                if i_obs is None:
                    corrIntens = obs.getLine(label=full_label).corrIntens
                else:
                    corrIntens = obs.getLine(label=full_label).corrIntens[i_obs]
                return corrIntens
            #pn.log_.debug('to eval is {0}'.format(eval_tem), calling=calling + 'TEST')
            try:
                value_tem = eval(eval_tem)
                #pn.log_.debug('to eval = {0}'.format(value_tem), calling=calling + 'TEST')
            except:
                pn.log_.warn('No value for {0} {1}: {2} in obs'.format(elem_tem, spec_tem, diag_tem), calling=calling)
                return None
        else:
            if type(value_tem) == type([]): value_tem = np.asarray(value_tem)
        if value_den is None:
            def L(wave): 
                if i_obs is None:
                    return obs.getLine(elem_den, spec_den, wave).corrIntens
                else:
                    return obs.getLine(elem_den, spec_den, wave).corrIntens[i_obs]
            def I(i, j):
                wave = atom_den.wave_Ang[i - 1, j - 1]
                pn.log_.debug('wave is {0}'.format(wave), calling=calling + 'TEST3')
                if i_obs is None:
                    return obs.getLine(elem_den, spec_den, wave).corrIntens
                else:
                    return obs.getLine(elem_den, spec_den, wave).corrIntens[i_obs]
            def B(label, I=I, L=L):
                full_label = elem_den + spec_den + '_' + label
                if i_obs is None:
                    corrIntens = obs.getLine(label=full_label).corrIntens
                else:
                    corrIntens = obs.getLine(label=full_label).corrIntens[i_obs]
                return corrIntens
            #pn.log_.debug('to eval is {0}'.format(eval_den), calling=calling + ' TEST')
            try:
                value_den = eval(eval_den)
                #pn.log_.debug('to eval = {0}'.format(value_den), calling=calling + ' TEST1')
            except:
                pn.log_.warn('No value for {0} {1}: {2} in obs'.format(elem_den, spec_den, diag_den), calling=calling)
                return None
        else:
            if type(value_den) == type([]): value_den = np.asarray(value_den)
        den = atom_den.getTemDen(value_den, tem=guess_tem, to_eval=eval_den,
                                 maxError=maxError, start_x=start_den, end_x=end_den)
        
        tem = atom_tem.getTemDen(value_tem, den=den, to_eval=eval_tem,
                                 maxError=maxError, start_x=start_tem, end_x=end_tem)
#        self.log_.debug('tem: ' + str(tem) + ' den:' + str(den), calling='getCrossTemDen')
        no_conv = np.ones_like(den).astype(bool)
        n_tot = np.asarray(value_tem).size
        for i in np.arange(max_iter):
            if type(tem) == type(1.):
                tem_old = tem
            else:
                tem_old = tem.copy()
            if type(den) == type(1.):
                den_old = den
            else:
                den_old = den.copy()
            
            if n_tot > 1:
                den[no_conv] = atom_den.getTemDen(value_den[no_conv], tem=tem_old[no_conv],
                                                  to_eval=eval_den, start_x=start_den, end_x=end_den)
                tem[no_conv] = atom_tem.getTemDen(value_tem[no_conv], den=den_old[no_conv],
                                                  to_eval=eval_tem, start_x=start_tem, end_x=end_tem)
            else:
                den = atom_den.getTemDen(value_den, tem=tem_old, to_eval=eval_den, start_x=start_den, end_x=end_den)
                tem = atom_tem.getTemDen(value_tem, den=den_old, to_eval=eval_tem, start_x=start_tem, end_x=end_tem)
                
            no_conv = ((abs(den_old - den) / den * 100) > tol_den) | ((abs(tem_old - tem) / tem * 100) > tol_tem)
            if type(no_conv) == type(True):
                n_no_conv = int(no_conv)
            else:
                n_no_conv = no_conv.sum()
            
            pn.log_.message('{0} (max={1}): not converged {2} of {3}.'.format(i, max_iter, n_no_conv, n_tot),
                            calling=calling)
            if n_no_conv == 0:
                return tem, den
        if n_tot == 1:
            tem = np.nan
            den = np.nan
        else:
            tem[no_conv] = np.nan
            den[no_conv] = np.nan
        return tem, den
    
    def getDiagLimits(self, diag):
        """
        Return the low and high density values for a given diagnostic
        """
        atom = self.atomDict[self.diags[diag][0]]
        to_eval = self.diags[diag][1]
        HDR = atom.getHighDensRatio(to_eval = to_eval)
        LDR = atom.getLowDensRatio(to_eval = to_eval)
        return(np.sort((LDR, HDR)))
    
    