# Data at [lat=  -4.98879105856447      , lon=   28.0026147968828      , k=          52 , time step=        1400 ] 
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
input.moment=[[         0],]
input.ncol=[[         4],]
input.pcols=[[         4],]
input.nver=[[        72],]
input.radius_max=[[  0.50000000000000002E-004],]
input.tair=[[  0.25143236121305881E+003,  0.25383982688822439E+003,  0.26105094678227255E+003,  0.26439278100322070E+003,  0.26850325516631881E+003,  0.26840733639492652E+003,  0.26571074356296373E+003,  0.26168781957191635E+003,  0.25462663783799380E+003,  0.24615682495556598E+003,  0.23800180507769215E+003,  0.23431359900233798E+003,  0.23188662948394759E+003,  0.22898120391932090E+003,  0.22296945564416848E+003,  0.21579655814937499E+003,  0.20968658086499204E+003,  0.20555944399763769E+003,  0.20126316692897527E+003,  0.19749133957256734E+003,  0.19357285613706910E+003,  0.19091419748166172E+003,  0.18983696817384543E+003,  0.18999641866883405E+003,  0.19047841336952459E+003,  0.19214948854555453E+003,  0.19112667192999152E+003,  0.18891065268113616E+003,  0.19132521374302948E+003,  0.19534080592527243E+003,  0.19915290105468048E+003,  0.20294970067986188E+003,  0.20733082619987269E+003,  0.21198693198028548E+003,  0.21660308866949219E+003,  0.22114213119213809E+003,  0.22569590699079296E+003,  0.23019280773565262E+003,  0.23480548868186213E+003,  0.23965054422568670E+003,  0.24429769341959403E+003,  0.24871376974252652E+003,  0.25356937869886573E+003,  0.25767261309325158E+003,  0.26100009366371364E+003,  0.26400283744806654E+003,  0.26670456736492446E+003,  0.26868731076659441E+003,  0.27025017176030792E+003,  0.27300215761107125E+003,  0.27525377490214362E+003,  0.27770764369092495E+003,  0.28108466479044904E+003,  0.28396718568415770E+003,  0.28633718041295685E+003,  0.28790789322254250E+003,  0.28905308482275160E+003,  0.29001124833810638E+003,  0.29084976345296536E+003,  0.29160029186182561E+003,  0.29227885493216729E+003,  0.29286808543617065E+003,  0.29338513973510095E+003,  0.29388117882770939E+003,  0.29423542939321248E+003,  0.29431618687679651E+003,  0.29378471992852729E+003,  0.29347084783559859E+003,  0.29408063515400812E+003,  0.29481401687193875E+003,  0.29552701271291829E+003,  0.29609528017143850E+003,],]
input.pmid=[[  0.12382541305561677E+002,  0.18282923550684863E+002,  0.26994886212093746E+002,  0.39858170362288227E+002,  0.58850914656475304E+002,  0.86893857004050020E+002,  0.12829949082549251E+003,  0.18943524794063876E+003,  0.27970269352932593E+003,  0.41298331550248361E+003,  0.59684493671180769E+003,  0.83774043789417703E+003,  0.11473787230534558E+004,  0.15333938222216857E+004,  0.19996337978167076E+004,  0.25444696509961595E+004,  0.31593251291285251E+004,  0.38366283094427008E+004,  0.45671197939426802E+004,  0.53309561434392535E+004,  0.61015181670369457E+004,  0.68476390229712551E+004,  0.75355335897260984E+004,  0.81946275124734711E+004,  0.88910543149560908E+004,  0.96466673444937041E+004,  0.10466496728675462E+005,  0.11356000142146328E+005,  0.12321099142746600E+005,  0.13368217993210610E+005,  0.14504326553270042E+005,  0.15736987784007928E+005,  0.17074407952158203E+005,  0.18474678442333348E+005,  0.19922184022484762E+005,  0.21471750653618299E+005,  0.23153008404824875E+005,  0.24977149267213306E+005,  0.26956316401515789E+005,  0.29103684667262849E+005,  0.31433548211621535E+005,  0.33961416456917555E+005,  0.36704117880872058E+005,  0.39679910706224066E+005,  0.42908603785385538E+005,  0.46362302456389851E+005,  0.49915621940473779E+005,  0.53452106382349368E+005,  0.56900916605652892E+005,  0.60227928346052591E+005,  0.63419041272566697E+005,  0.66431314882017032E+005,  0.69222574753536552E+005,  0.71752382839743819E+005,  0.73983010292934385E+005,  0.75880372454932469E+005,  0.77443028053749105E+005,  0.78764362379486440E+005,  0.79990366197702984E+005,  0.81187255281842212E+005,  0.82353108205777855E+005,  0.83486033117299841E+005,  0.84584176201802504E+005,  0.85645723135198685E+005,  0.86668898758958749E+005,  0.87651971686794801E+005,  0.88593266031427993E+005,  0.89491172888908550E+005,  0.90344146839537454E+005,  0.91155338238683689E+005,  0.91850554035594876E+005,  0.92288613013553535E+005,],]
input.radius_part=[[  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.15109629247903635E-005,  0.19407136234505346E-005,  0.22442677634572075E-005,  0.25525179549420406E-005,  0.28417343287006439E-005,  0.30717728411024080E-005,  0.32187079430850649E-005,  0.33295108568805802E-005,  0.33581389134597751E-005,  0.26095991798367073E-005,  0.14746357036944866E-005,  0.10027790966036392E-005,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.86563926892070719E-006,  0.89935209802971794E-006,  0.90559557351473757E-006,  0.91774689366844418E-006,  0.94131548879414300E-006,  0.95791677900665382E-006,  0.96469363885711524E-006,  0.96341473931759348E-006,  0.96078604540569116E-006,  0.95186602522176238E-006,  0.93632102899969262E-006,  0.92685737319478648E-006,  0.91865638203538930E-006,  0.91670344820442098E-006,  0.92310410605960610E-006,  0.93187684505460990E-006,  0.10997883418000315E-005,  0.96364897053197941E-006,  0.95274110775058209E-006,  0.94535331831377239E-006,  0.20971660065521869E-005,  0.17593006354496512E-005,  0.25086241731402559E-005,  0.27789800957680936E-005,  0.17346479358692447E-005,  0.13269612565733259E-005,  0.11713697437609708E-005,  0.11041648399691404E-005,  0.10862959814520149E-005,  0.10822016321083548E-005,  0.10832956110022361E-005,  0.10877025772770995E-005,  0.10937292709406852E-005,  0.11016680008491996E-005,  0.11146061463783655E-005,  0.11281259698468912E-005,  0.11420292638804319E-005,  0.11505451032884557E-005,  0.11728179328135246E-005,  0.13430247455527368E-005,  0.17998345629096252E-005,  0.18004692444715021E-005,  0.16900697747443475E-005,  0.16012641136138683E-005,  0.15530448837092021E-005,],]
input.density_part=[[  0.10000000000000000E+001,  0.17699904624872529E+004,  0.17699999659878092E+004,  0.17699999997457764E+004,  0.17699999999995503E+004,  0.17699999999999993E+004,  0.17700000000000000E+004,  0.17700000000000000E+004,  0.17700000000000000E+004,  0.17700000000000000E+004,  0.17700000000000002E+004,  0.17700000000001271E+004,  0.17700000000261459E+004,  0.17700000026598220E+004,  0.17700001661304641E+004,  0.17700174949055704E+004,  0.17711035762570175E+004,  0.17849359207653822E+004,  0.18444295626719068E+004,  0.18704330054275933E+004,  0.18777770282395390E+004,  0.18708305348690562E+004,  0.17638846465805329E+004,  0.18090356623985672E+004,  0.19104632314448211E+004,  0.19194047692958773E+004,  0.18483860325811115E+004,  0.17705271262302335E+004,  0.17760411002102312E+004,  0.17933871063340464E+004,  0.18297388342728693E+004,  0.18492711579558961E+004,  0.18535434137419024E+004,  0.18566092491933493E+004,  0.18651459808462812E+004,  0.18695885361894505E+004,  0.18783090619188401E+004,  0.18925786677088556E+004,  0.19124816878029094E+004,  0.19303573261055594E+004,  0.19333894232828236E+004,  0.19120291336758162E+004,  0.18902639737029112E+004,  0.15432685741501787E+004,  0.18105805473632010E+004,  0.18249391188392210E+004,  0.18349088073633366E+004,  0.10727914979764500E+004,  0.11246261777782568E+004,  0.10541932899919329E+004,  0.10472446061501967E+004,  0.12127484730448286E+004,  0.14738992055433630E+004,  0.16838682808556139E+004,  0.18048020988071810E+004,  0.18176765620756753E+004,  0.18059398269679152E+004,  0.17921332237334705E+004,  0.17779374014420571E+004,  0.17612806444055043E+004,  0.17402660952886467E+004,  0.17060208920850694E+004,  0.16711481244933771E+004,  0.16424792967494604E+004,  0.16239102287810560E+004,  0.15798498517030359E+004,  0.13510018555247950E+004,  0.11369007874785329E+004,  0.11343042654553801E+004,  0.11609070962071123E+004,  0.11882312517683481E+004,  0.12058645526383575E+004,],]
input.sig_part=[[  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,],]
# Output data is stored here.
output = Object()
output.vlc_grv=[[  0.98010376860593003E-004,  0.11807374244524423E+000,  0.81115929007633772E-001,  0.55308268184911726E-001,  0.37768575769354135E-001,  0.25595717790279117E-001,  0.31174224661828506E-001,  0.27084082293660609E-001,  0.21139662897730838E-001,  0.16302913182760247E-001,  0.12716289219452585E-001,  0.10118876299880683E-001,  0.81194717638835810E-002,  0.66902842104937955E-002,  0.55612074643388403E-002,  0.33767678744048750E-002,  0.13862900918504217E-002,  0.74642683148722164E-003,  0.53966009175245933E-003,  0.48442196689062487E-003,  0.43984758806507731E-003,  0.40447339709476096E-003,  0.35821103961479608E-003,  0.34861498877827217E-003,  0.35051859521805379E-003,  0.33571289553634175E-003,  0.30926084922726604E-003,  0.29779287251376746E-003,  0.30358471642526094E-003,  0.29700863498524638E-003,  0.29705444191245020E-003,  0.30051470106171258E-003,  0.29800969837260211E-003,  0.29094796947147771E-003,  0.28164302943925842E-003,  0.27170372093490735E-003,  0.25997785765548108E-003,  0.24664996950404030E-003,  0.23723766782567676E-003,  0.22843264734426338E-003,  0.22100094011752214E-003,  0.21469799335624898E-003,  0.20943772354998968E-003,  0.22300127991659819E-003,  0.20191843719610439E-003,  0.19485969122491812E-003,  0.18916263841929948E-003,  0.47631810477317374E-003,  0.35425347988884222E-003,  0.64513352188562943E-003,  0.77317410456903828E-003,  0.35928649774072602E-003,  0.26011535468827078E-003,  0.23279811107240129E-003,  0.22145874945553640E-003,  0.21488262846738337E-003,  0.21087969486549807E-003,  0.20874217163164408E-003,  0.20784652380667696E-003,  0.20728103619745803E-003,  0.20689634474275068E-003,  0.20666659216240141E-003,  0.20648392385245323E-003,  0.20710540306037844E-003,  0.20717352485492545E-003,  0.20861990033386611E-003,  0.23038703471830492E-003,  0.33885957350632202E-003,  0.33757720818982647E-003,  0.30525476503996614E-003,  0.28112226708946309E-003,  0.26863772991593050E-003,],]