# Data at [lat=  -4.98879105856447      , lon=   28.0026147968828      , k=          52 , time step=        1388 ] 
# This file was generated by E3SM.

from math import nan as nan, inf as inf

# Object is just a dynamic container that stores input/output data.
class Object(object):
    pass
# Settings are stored here.
settings = Object()
# Input is stored here.
input = Object()
input.dt = [    0.00, ]
input.moment=[[         3],]
input.ncol=[[         4],]
input.pcols=[[         4],]
input.nver=[[        72],]
input.radius_max=[[  0.50000000000000002E-004],]
input.tair=[[  0.25318555702947307E+003,  0.26346079788871816E+003,  0.26529643429967541E+003,  0.26292479014276523E+003,  0.26492640909518741E+003,  0.26713115839323200E+003,  0.26485881969255104E+003,  0.25900640611258183E+003,  0.24845116962080573E+003,  0.24013513139192264E+003,  0.23594691863020864E+003,  0.23594636462882659E+003,  0.23463235118128634E+003,  0.23144822552072054E+003,  0.22409302469104605E+003,  0.21584898464029027E+003,  0.21006226543087715E+003,  0.20621513088826612E+003,  0.20163419510415102E+003,  0.19707692185778032E+003,  0.19217872352382344E+003,  0.18919155308225154E+003,  0.18824305259886270E+003,  0.18830309539303519E+003,  0.18856745625607249E+003,  0.18882372743518872E+003,  0.18677022669794331E+003,  0.18802740712788696E+003,  0.19133491326552118E+003,  0.19508481040438991E+003,  0.19931069308333494E+003,  0.20356419510073565E+003,  0.20797316901158115E+003,  0.21234598132945283E+003,  0.21671629644385445E+003,  0.22123631674532294E+003,  0.22591525220990704E+003,  0.23056996296996400E+003,  0.23522079389933967E+003,  0.23981163827495016E+003,  0.24436980846980933E+003,  0.24832677064235551E+003,  0.25248433644059165E+003,  0.25699120611339583E+003,  0.26181887248847448E+003,  0.26560568647660449E+003,  0.26788029057956095E+003,  0.27029536335129626E+003,  0.27139606423566437E+003,  0.27323284748095642E+003,  0.27539166126815752E+003,  0.27896232078256816E+003,  0.28163368925804457E+003,  0.28400938856411244E+003,  0.28605002823956028E+003,  0.28762899756631799E+003,  0.28884375073724908E+003,  0.28980636897875530E+003,  0.29064696927919255E+003,  0.29143228080267903E+003,  0.29221728640491938E+003,  0.29291127763945178E+003,  0.29350138516001482E+003,  0.29401248929861526E+003,  0.29447215408416588E+003,  0.29492384586096699E+003,  0.29543429129897328E+003,  0.29595523959982518E+003,  0.29644876242257982E+003,  0.29494377899669843E+003,  0.29386300564920026E+003,  0.29372011965820894E+003,],]
input.pmid=[[  0.12382541305561677E+002,  0.18282923550684863E+002,  0.26994886212093746E+002,  0.39858170362288227E+002,  0.58850914656475304E+002,  0.86893857004050020E+002,  0.12829949082549251E+003,  0.18943524794063876E+003,  0.27970269352932593E+003,  0.41298331550248361E+003,  0.59684493671180769E+003,  0.83774043789417703E+003,  0.11473787230534558E+004,  0.15333938222216857E+004,  0.19996337978167076E+004,  0.25444696509961595E+004,  0.31593251291285251E+004,  0.38366283094427008E+004,  0.45671197939426802E+004,  0.53309561434392535E+004,  0.61015181670369457E+004,  0.68476390229712551E+004,  0.75355335897260984E+004,  0.81946275124734711E+004,  0.88910543149560908E+004,  0.96466673444937041E+004,  0.10466496728675462E+005,  0.11356000142146328E+005,  0.12321099142746600E+005,  0.13368217993210610E+005,  0.14504326553270042E+005,  0.15736987784007928E+005,  0.17074407952158203E+005,  0.18474134462859154E+005,  0.19920281499644407E+005,  0.21468149776477287E+005,  0.23147564836852936E+005,  0.24969706405782075E+005,  0.26946704334881644E+005,  0.29091719044027304E+005,  0.31419029014330416E+005,  0.33944126668413403E+005,  0.36683822039046914E+005,  0.39656353338680070E+005,  0.42881507709325233E+005,  0.46331421061076311E+005,  0.49880846038598727E+005,  0.53413454425494238E+005,  0.56858484686969066E+005,  0.60181849959001971E+005,  0.63369465365441021E+005,  0.66378437466228163E+005,  0.69166638063830789E+005,  0.71693673431735500E+005,  0.73921856074347234E+005,  0.75817138690852415E+005,  0.77378081589019392E+005,  0.78697967707180782E+005,  0.79922627802553063E+005,  0.81118205074029494E+005,  0.82282780201103407E+005,  0.83414463405405564E+005,  0.84511402903589886E+005,  0.85571786360807746E+005,  0.86593840564906772E+005,  0.87575836026372548E+005,  0.88516098693964625E+005,  0.89413021427973537E+005,  0.90265060503291883E+005,  0.91075417542213938E+005,  0.91769991565216711E+005,  0.92207666319979267E+005,],]
input.radius_part=[[  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.14878566713285527E-005,  0.19340574296112383E-005,  0.23065539766526600E-005,  0.26436863185593783E-005,  0.28855301421480876E-005,  0.30564698012937662E-005,  0.31856897971361915E-005,  0.33118889037508122E-005,  0.33581389134597751E-005,  0.26981865339361036E-005,  0.14536624237081902E-005,  0.10031335944457302E-005,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.84068512593764513E-006,  0.86103875377950829E-006,  0.84233564496515528E-006,  0.83953472836494377E-006,  0.84923232353919449E-006,  0.84416927220134823E-006,  0.87924932680672093E-006,  0.89256156919785289E-006,  0.91567459797568703E-006,  0.93883528401297040E-006,  0.94668543901646821E-006,  0.95299223313648000E-006,  0.95916433108092506E-006,  0.96257566400362882E-006,  0.96459236477208693E-006,  0.96257309959255253E-006,  0.95828846795602744E-006,  0.95679234070977724E-006,  0.95208283119633592E-006,  0.95003405575494531E-006,  0.94628977183990792E-006,  0.94471490981789270E-006,  0.94865986628187532E-006,  0.12536484281553406E-005,  0.13272384251694306E-005,  0.14844349863236396E-005,  0.14934222070332987E-005,  0.15706790605379101E-005,  0.19928268348929103E-005,  0.19477875812576890E-005,  0.14218384273525824E-005,  0.12516349901658104E-005,  0.11567168259033397E-005,  0.11019118401637568E-005,  0.10896747828004306E-005,  0.10913664333583550E-005,  0.10917560864000799E-005,  0.10942956412247215E-005,  0.10971869698710206E-005,  0.10977429462116007E-005,  0.11018342793552813E-005,  0.11111224135797487E-005,  0.11285797625811941E-005,  0.11626261944518993E-005,  0.12306549425696091E-005,  0.12913970850697840E-005,  0.13225853908506714E-005,  0.13135743732642630E-005,  0.15577953876808743E-005,  0.20801795487469263E-005,  0.23753671399326282E-005,],]
input.density_part=[[  0.10000000000000000E+001,  0.17699979350898268E+004,  0.17699999896368790E+004,  0.17699999999747483E+004,  0.17699999999998593E+004,  0.17699999999999993E+004,  0.17700000000000000E+004,  0.17700000000000000E+004,  0.17700000000000000E+004,  0.17700000000000000E+004,  0.17700000000000002E+004,  0.17700000000000234E+004,  0.17700000000043749E+004,  0.17700000011988818E+004,  0.17700002783115722E+004,  0.17700337361962595E+004,  0.17715094581419523E+004,  0.17865586267899637E+004,  0.18413342309584384E+004,  0.18674669274804428E+004,  0.18828222605953868E+004,  0.18501542708077657E+004,  0.17309009775823245E+004,  0.18729066973388522E+004,  0.19305295377709649E+004,  0.18987742266401772E+004,  0.18230539561267071E+004,  0.17667479821047596E+004,  0.17915052550397165E+004,  0.18296832532920055E+004,  0.18397934593936873E+004,  0.18346021659239041E+004,  0.18311728678149886E+004,  0.18383341029743510E+004,  0.18485132628562342E+004,  0.18554266130654373E+004,  0.18610850206554505E+004,  0.18687987441263529E+004,  0.18732773914900022E+004,  0.18841765065281543E+004,  0.18901180946948964E+004,  0.18911171497313239E+004,  0.19039137766910596E+004,  0.19036174351254817E+004,  0.13959953837437097E+004,  0.13420626419249827E+004,  0.12285323436435715E+004,  0.12120438219692744E+004,  0.12012896409411712E+004,  0.11180834129829514E+004,  0.11388327746418436E+004,  0.13622119139846966E+004,  0.15292119202726935E+004,  0.16627541078219253E+004,  0.17469756830846134E+004,  0.17504416121524735E+004,  0.17412607473047283E+004,  0.17420981736990391E+004,  0.17431730942344038E+004,  0.17446514921632574E+004,  0.17510681781605056E+004,  0.17434018653365147E+004,  0.17240042446540267E+004,  0.16908362460820501E+004,  0.16269359485126104E+004,  0.15152225248520213E+004,  0.14374950156570110E+004,  0.14092821359430579E+004,  0.14218271728489087E+004,  0.12414353488196909E+004,  0.10941429952255876E+004,  0.10625201474169605E+004,],]
input.sig_part=[[  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,],]
# Output data is stored here.
output = Object()
output.vlc_grv=[[  0.27747223011288265E-003,  0.33945589406583837E+000,  0.23086639375129686E+000,  0.15583019587871466E+000,  0.10610352018680394E+000,  0.72321729861719267E-001,  0.87466367312121443E-001,  0.77603679341547827E-001,  0.63500806145339428E-001,  0.51346784847207838E-001,  0.41763264678528957E-001,  0.34998657144859051E-001,  0.30363538674336488E-001,  0.27711111247393897E-001,  0.25614679128754998E-001,  0.16636945338844425E-001,  0.57226158630376200E-002,  0.30026539906139781E-002,  0.21807984640368136E-002,  0.20695681085148901E-002,  0.19902554066824054E-002,  0.18908609338947713E-002,  0.17882927429986068E-002,  0.18171495816731820E-002,  0.18180988163798185E-002,  0.17834874624910193E-002,  0.16704830431650056E-002,  0.17013707583980024E-002,  0.17239968354577652E-002,  0.17928779136945319E-002,  0.18334345446668324E-002,  0.18056594806266867E-002,  0.17752417055099980E-002,  0.17573054439941648E-002,  0.17352759447164111E-002,  0.17063669303965454E-002,  0.16639476035322793E-002,  0.16181947795627424E-002,  0.15803582102177068E-002,  0.15399925296624008E-002,  0.15054337887977895E-002,  0.14659985761785186E-002,  0.14426011312152504E-002,  0.14248874508433294E-002,  0.17465262744173339E-002,  0.18451619727526533E-002,  0.20760754392882366E-002,  0.20510721360254697E-002,  0.22291699272415145E-002,  0.32822986836522996E-002,  0.31715359791798818E-002,  0.20242729888838339E-002,  0.17560884971128468E-002,  0.16245464041698932E-002,  0.15424249885207435E-002,  0.15041644751359713E-002,  0.14945593489966646E-002,  0.14913558023507130E-002,  0.14946692740860543E-002,  0.14995352551240735E-002,  0.15024888013275026E-002,  0.15031497963520712E-002,  0.15077473519434219E-002,  0.15215268711467151E-002,  0.15488243998633972E-002,  0.16093909623194061E-002,  0.16748881803464659E-002,  0.17176451428606025E-002,  0.17071233834990132E-002,  0.20902021660722551E-002,  0.32652521944809094E-002,  0.41227096552944294E-002,],]