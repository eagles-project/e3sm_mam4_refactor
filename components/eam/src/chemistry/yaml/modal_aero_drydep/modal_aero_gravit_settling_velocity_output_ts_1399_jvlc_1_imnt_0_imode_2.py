# Data at [lat=  -4.98879105856447      , lon=   28.0026147968828      , k=          52 , time step=        1399 ] 
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
input.tair=[[  0.25106150734556880E+003,  0.25358590510502427E+003,  0.26088053377168751E+003,  0.26412591401951693E+003,  0.26831783648551072E+003,  0.26835232090573948E+003,  0.26540598605688638E+003,  0.26122190002652911E+003,  0.25415494810369276E+003,  0.24569465978019412E+003,  0.23755139621200831E+003,  0.23400649595550598E+003,  0.23162528504630330E+003,  0.22871336865182977E+003,  0.22266886288884729E+003,  0.21554317174431526E+003,  0.20962424792050180E+003,  0.20564997566190289E+003,  0.20139758995287994E+003,  0.19765142094376955E+003,  0.19370106176836683E+003,  0.19095115516058536E+003,  0.18977841092785249E+003,  0.18986605616574491E+003,  0.19029492680366619E+003,  0.19179174013810970E+003,  0.19033270867461476E+003,  0.18852605745373182E+003,  0.19133605694726464E+003,  0.19531796806056434E+003,  0.19913473679182295E+003,  0.20293390432945534E+003,  0.20728510997707045E+003,  0.21190952476290721E+003,  0.21651996722397774E+003,  0.22106342037868947E+003,  0.22564734231229332E+003,  0.23015000800985612E+003,  0.23468705368684505E+003,  0.23944999726186546E+003,  0.24412035207508578E+003,  0.24858553004218786E+003,  0.25345035043835335E+003,  0.25761175013406626E+003,  0.26094953734937769E+003,  0.26401029108225873E+003,  0.26678474522190453E+003,  0.26880889498660815E+003,  0.27021476055419873E+003,  0.27291232080672040E+003,  0.27520231164056753E+003,  0.27770227017222891E+003,  0.28106096539122592E+003,  0.28394827328696306E+003,  0.28631104361705013E+003,  0.28782587223269815E+003,  0.28899155644261214E+003,  0.28995345645037958E+003,  0.29078669200131134E+003,  0.29154074271743946E+003,  0.29222182287361710E+003,  0.29280943507700340E+003,  0.29334081808190547E+003,  0.29384527051102560E+003,  0.29419541315209966E+003,  0.29438193499561231E+003,  0.29430944346660215E+003,  0.29335755817498426E+003,  0.29279030584250899E+003,  0.29329781874223119E+003,  0.29389879217947407E+003,  0.29433696772563007E+003,],]
input.pmid=[[  0.12382541305561677E+002,  0.18282923550684863E+002,  0.26994886212093746E+002,  0.39858170362288227E+002,  0.58850914656475304E+002,  0.86893857004050020E+002,  0.12829949082549251E+003,  0.18943524794063876E+003,  0.27970269352932593E+003,  0.41298331550248361E+003,  0.59684493671180769E+003,  0.83774043789417703E+003,  0.11473787230534558E+004,  0.15333938222216857E+004,  0.19996337978167076E+004,  0.25444696509961595E+004,  0.31593251291285251E+004,  0.38366283094427008E+004,  0.45671197939426802E+004,  0.53309561434392535E+004,  0.61015181670369457E+004,  0.68476390229712551E+004,  0.75355335897260984E+004,  0.81946275124734711E+004,  0.88910543149560908E+004,  0.96466673444937041E+004,  0.10466496728675462E+005,  0.11356000142146328E+005,  0.12321099142746600E+005,  0.13368217993210610E+005,  0.14504326553270042E+005,  0.15736987784007928E+005,  0.17074407952158203E+005,  0.18474845077873411E+005,  0.19922766816387884E+005,  0.21472853699110929E+005,  0.23154675916067310E+005,  0.24979429215801516E+005,  0.26959260835985620E+005,  0.29107350059404700E+005,  0.31437995832159912E+005,  0.33966712784043179E+005,  0.36710335042959960E+005,  0.39687126961339134E+005,  0.42916904042071830E+005,  0.46371762258147231E+005,  0.49926274734310718E+005,  0.53463946516050310E+005,  0.56913914643301730E+005,  0.60242043395012413E+005,  0.63434227705970668E+005,  0.66447512656299135E+005,  0.69239709665680508E+005,  0.71770367110059742E+005,  0.74001743474469346E+005,  0.75899742657169001E+005,  0.77462922902247607E+005,  0.78784700852906448E+005,  0.80011116289561687E+005,  0.81208407217192987E+005,  0.82374651564666186E+005,  0.83507956844497996E+005,  0.84606468619853316E+005,  0.85668371957287745E+005,  0.86691891102043926E+005,  0.87675294086801310E+005,  0.88616904461853061E+005,  0.89515112782815166E+005,  0.90368373110889646E+005,  0.91179820097007745E+005,  0.91875232486523804E+005,  0.92313409162378841E+005,],]
input.radius_part=[[  0.13325053635595470E-007,  0.13990944496573431E-007,  0.17191877826176230E-007,  0.18851609635316096E-007,  0.19981491524779117E-007,  0.23927830904270331E-007,  0.27302291250127147E-007,  0.29326659082689562E-007,  0.31813720497647043E-007,  0.33920776648049540E-007,  0.35496487878461187E-007,  0.36340085112598469E-007,  0.37074498079950922E-007,  0.38026141416642339E-007,  0.37909106715773882E-007,  0.37304722648963576E-007,  0.38406478355327668E-007,  0.40018633454294873E-007,  0.39597181414625468E-007,  0.34848451530870955E-007,  0.30276112277324293E-007,  0.28035529788343221E-007,  0.27508929003912788E-007,  0.27453786754145945E-007,  0.26441973303017286E-007,  0.21785631172753751E-007,  0.14645021092567291E-007,  0.15081877599368940E-007,  0.17279788014208834E-007,  0.19050417692577246E-007,  0.20925740601953762E-007,  0.22302360584227038E-007,  0.25133162705951066E-007,  0.28989769481066932E-007,  0.30657988502529241E-007,  0.30982134638696878E-007,  0.31620016431971417E-007,  0.32389219274138733E-007,  0.32866586394352394E-007,  0.33011077170558169E-007,  0.33045971342942383E-007,  0.32473909823037626E-007,  0.32632045632805300E-007,  0.35578554816865148E-007,  0.34739195137337482E-007,  0.35446539090320303E-007,  0.38277166597857387E-007,  0.44562559974892381E-007,  0.48583957699266557E-007,  0.54539448070992928E-007,  0.56439917709658104E-007,  0.40963208659847442E-007,  0.33613381889853219E-007,  0.30508895369589904E-007,  0.28602332035635669E-007,  0.26166789121631027E-007,  0.23382025338203618E-007,  0.22012846381973297E-007,  0.20497080749385901E-007,  0.19338295895965236E-007,  0.18512407911177341E-007,  0.18072541610285444E-007,  0.17957553470342239E-007,  0.18033165444003090E-007,  0.18401797877866925E-007,  0.20383515067376445E-007,  0.23041190603960513E-007,  0.31249052474077969E-007,  0.46681113273038304E-007,  0.45971250112719421E-007,  0.41358134763516017E-007,  0.39208536318330984E-007,],]
input.density_part=[[  0.17711753672442107E+004,  0.17717897128906388E+004,  0.17718408080468698E+004,  0.17718833638353024E+004,  0.17718600558730882E+004,  0.17718696991122717E+004,  0.17714300616927801E+004,  0.17128489929015088E+004,  0.16095177213353813E+004,  0.15696382031309342E+004,  0.15512357878816420E+004,  0.15450598832330763E+004,  0.15389841096325149E+004,  0.15343351245439408E+004,  0.15181414746428118E+004,  0.14532256723035434E+004,  0.13400077291890318E+004,  0.12591502125112120E+004,  0.12247748954032866E+004,  0.11963306033993533E+004,  0.11836488194916626E+004,  0.11737525208974425E+004,  0.11567134127462243E+004,  0.11499983784961105E+004,  0.11505497637383351E+004,  0.11599662021044521E+004,  0.12186924143504052E+004,  0.11766361019363981E+004,  0.11572522802030871E+004,  0.11352946468441380E+004,  0.11084992630944207E+004,  0.10699351624988285E+004,  0.10351376772581623E+004,  0.10215605019157898E+004,  0.10197720401568344E+004,  0.10202442660036115E+004,  0.10203580767507965E+004,  0.10200942753200213E+004,  0.10220803678588362E+004,  0.10360884395821195E+004,  0.10527908540356016E+004,  0.10615047042502281E+004,  0.10626647033957611E+004,  0.10521818474911047E+004,  0.10636957362603748E+004,  0.10618612493945548E+004,  0.10489421761363290E+004,  0.10356639966331763E+004,  0.10396681451735169E+004,  0.10376435461629630E+004,  0.10425464038663144E+004,  0.11265416479099970E+004,  0.12336080139714129E+004,  0.13118334642576617E+004,  0.13611104657777066E+004,  0.13948692427024371E+004,  0.14211846605088544E+004,  0.14436601105758857E+004,  0.14861187494868482E+004,  0.15092748213149898E+004,  0.15095648896149228E+004,  0.14866931597530518E+004,  0.14839621940990341E+004,  0.14849100986187577E+004,  0.14930066396618602E+004,  0.14834434006116019E+004,  0.14049632537846003E+004,  0.12163612947206045E+004,  0.10638012704449466E+004,  0.10665102057922134E+004,  0.10911832447140841E+004,  0.11069478153524615E+004,],]
input.sig_part=[[  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,  0.16000000000000001E+001,],]
# Output data is stored here.
output = Object()
output.vlc_grv=[[  0.25859758240008206E-002,  0.18488065733518674E-002,  0.15606531270203499E-002,  0.11662612637159666E-002,  0.84384040038473860E-003,  0.68445717875209448E-003,  0.52592938775197166E-003,  0.36706155749216546E-003,  0.25000231640091798E-003,  0.17314992560479923E-003,  0.12188094156784937E-003,  0.87922227391464656E-004,  0.64942771771510648E-004,  0.49419137229663555E-004,  0.36922047534934109E-004,  0.26925050409326853E-004,  0.20336620960337518E-004,  0.16274462611563430E-004,  0.13046685967395035E-004,  0.95243572589024361E-005,  0.70830698237930629E-005,  0.57582224217894052E-005,  0.50499110984421613E-005,  0.46140108336116315E-005,  0.41050400466433820E-005,  0.31493618960416322E-005,  0.20348442811376788E-005,  0.18583765281913810E-005,  0.19493222130984071E-005,  0.19677721954830181E-005,  0.19692574398359591E-005,  0.18894189560682651E-005,  0.19266981464973006E-005,  0.20606499906425152E-005,  0.20460549224312058E-005,  0.19433505122667603E-005,  0.18633467895448616E-005,  0.17917706443204982E-005,  0.17093985015464102E-005,  0.16324522074182602E-005,  0.15564047056652727E-005,  0.14430104835765518E-005,  0.13604518718770945E-005,  0.13815980882414830E-005,  0.12724609877790028E-005,  0.12136768216709617E-005,  0.12221013644781789E-005,  0.13437012503150296E-005,  0.14079800196303494E-005,  0.15303380927108952E-005,  0.15348019152178602E-005,  0.11100109829141561E-005,  0.94526686995239396E-006,  0.87851346874360660E-006,  0.82892882158977650E-006,  0.75525135394358023E-006,  0.67007409540678185E-006,  0.62897443815237599E-006,  0.59212558236249755E-006,  0.55805728396985595E-006,  0.52633667568376141E-006,  0.49929374374480549E-006,  0.48934523083975014E-006,  0.48646236857523521E-006,  0.49443938950484501E-006,  0.54226814070660465E-006,  0.58009706246758452E-006,  0.69405015273462815E-006,  0.95087694017467129E-006,  0.92987989002284806E-006,  0.83653924136614418E-006,  0.79524605941861843E-006,],]