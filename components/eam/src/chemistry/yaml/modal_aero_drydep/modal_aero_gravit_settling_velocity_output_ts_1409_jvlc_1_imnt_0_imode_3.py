# Data at [lat=  -4.98879105856447      , lon=   28.0026147968828      , k=          52 , time step=        1409 ] 
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
input.tair=[[  0.25448531923775388E+003,  0.26917507711550218E+003,  0.27527995484790199E+003,  0.27646320368299450E+003,  0.27172252748854567E+003,  0.26670929531280802E+003,  0.26359229089452674E+003,  0.25827808388259535E+003,  0.24999793370358381E+003,  0.24194112389804837E+003,  0.23612810366470794E+003,  0.23404715009753059E+003,  0.23239424896390403E+003,  0.23000030395266378E+003,  0.22488400259528524E+003,  0.21831966249088592E+003,  0.21125272062882010E+003,  0.20500071905936795E+003,  0.19934742129278064E+003,  0.19459643082999108E+003,  0.18991880780124160E+003,  0.18745431489477437E+003,  0.18720657966674224E+003,  0.18808602612324276E+003,  0.18955005485861017E+003,  0.19189372140924138E+003,  0.19246292185821773E+003,  0.18942419619735230E+003,  0.19092052060648689E+003,  0.19530521835040020E+003,  0.19895372900264405E+003,  0.20312470143001036E+003,  0.20763762754290809E+003,  0.21223112923394919E+003,  0.21679596350125865E+003,  0.22145864194805193E+003,  0.22621343845600012E+003,  0.23091885265942648E+003,  0.23553714421331264E+003,  0.24008905161479922E+003,  0.24440276511803674E+003,  0.24854858448827110E+003,  0.25276771067008607E+003,  0.25693665711977485E+003,  0.26103863891854525E+003,  0.26516762471027391E+003,  0.26765853652637441E+003,  0.26962471703376048E+003,  0.27146797654160008E+003,  0.27352337998539548E+003,  0.27532258695030811E+003,  0.27878137396870255E+003,  0.28147124890647115E+003,  0.28392761488813881E+003,  0.28610158929374558E+003,  0.28792231667936193E+003,  0.28932691344681268E+003,  0.29039396193302127E+003,  0.29128334984127679E+003,  0.29213778156204290E+003,  0.29291487406086998E+003,  0.29361045287862811E+003,  0.29419459517035313E+003,  0.29470365594657596E+003,  0.29494991964494443E+003,  0.29517145660896813E+003,  0.29563598263525000E+003,  0.29621059631874215E+003,  0.29674700642431537E+003,  0.29684167922222639E+003,  0.29552990453632617E+003,  0.29443385975310053E+003,],]
input.pmid=[[  0.12382541305561677E+002,  0.18282923550684863E+002,  0.26994886212093746E+002,  0.39858170362288227E+002,  0.58850914656475304E+002,  0.86893857004050020E+002,  0.12829949082549251E+003,  0.18943524794063876E+003,  0.27970269352932593E+003,  0.41298331550248361E+003,  0.59684493671180769E+003,  0.83774043789417703E+003,  0.11473787230534558E+004,  0.15333938222216857E+004,  0.19996337978167076E+004,  0.25444696509961595E+004,  0.31593251291285251E+004,  0.38366283094427008E+004,  0.45671197939426802E+004,  0.53309561434392535E+004,  0.61015181670369457E+004,  0.68476390229712551E+004,  0.75355335897260984E+004,  0.81946275124734711E+004,  0.88910543149560908E+004,  0.96466673444937041E+004,  0.10466496728675462E+005,  0.11356000142146328E+005,  0.12321099142746600E+005,  0.13368217993210610E+005,  0.14504326553270042E+005,  0.15736987784007928E+005,  0.17074407952158203E+005,  0.18473058038959040E+005,  0.19916516796832941E+005,  0.21461024378553546E+005,  0.23136793131402825E+005,  0.24954978508297114E+005,  0.26927684024841717E+005,  0.29068041527794463E+005,  0.31390298498353026E+005,  0.33909913720210810E+005,  0.36643660727361421E+005,  0.39609738135991014E+005,  0.42827890127680330E+005,  0.46270313108365823E+005,  0.49812031655515690E+005,  0.53336970123936953E+005,  0.56774520615591464E+005,  0.60090670274007964E+005,  0.63271364804975536E+005,  0.66273803896706071E+005,  0.69055950818173136E+005,  0.71577499544778417E+005,  0.73800844408276869E+005,  0.75692012030418933E+005,  0.77249565844989003E+005,  0.78566586257143019E+005,  0.79788587400379591E+005,  0.80981568863429944E+005,  0.82143615492377197E+005,  0.83272841612557997E+005,  0.84367399460072804E+005,  0.85425480636368346E+005,  0.86445315780774035E+005,  0.87425179157101345E+005,  0.88363400347253671E+005,  0.89258375702591526E+005,  0.90108564851256990E+005,  0.90917270862283287E+005,  0.91610574946241948E+005,  0.92047489402039384E+005,],]
input.radius_part=[[  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.15584069658707119E-005,  0.20069051081717733E-005,  0.23181701733847726E-005,  0.26020939574678828E-005,  0.28606941852710954E-005,  0.30756804375378203E-005,  0.32258161157320836E-005,  0.33381525976067360E-005,  0.33581389134597751E-005,  0.26570499757529901E-005,  0.14902880995604765E-005,  0.97758197262443375E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.90423504963877321E-006,  0.90712887791772208E-006,  0.85701792405774681E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.84615077280310167E-006,  0.91047283792253698E-006,  0.91614521628542164E-006,  0.93319490013851746E-006,  0.93954455845696562E-006,  0.93377828261734880E-006,  0.92703268815966430E-006,  0.92368394452051957E-006,  0.92166785069401513E-006,  0.91900177259091229E-006,  0.91666425793928372E-006,  0.91981758444145823E-006,  0.92385448556235938E-006,  0.92994120374001515E-006,  0.93212712431657938E-006,  0.93449115864991273E-006,  0.10094017163865631E-005,  0.11873506677492249E-005,  0.12458705327619589E-005,  0.14760253039131313E-005,  0.15351393131083958E-005,  0.13735356849638722E-005,  0.15578670835783559E-005,  0.17574087429953275E-005,  0.14264744130375594E-005,  0.13103406062482110E-005,  0.12271417814105204E-005,  0.11644680506188408E-005,  0.11180409711374074E-005,  0.10896991989516353E-005,  0.10698323457311796E-005,  0.10602684528145979E-005,  0.10552970370889340E-005,  0.10562985008154571E-005,  0.10604962965449227E-005,  0.10611466580100462E-005,  0.10699197633724115E-005,  0.11225875738966767E-005,  0.12144289640826848E-005,  0.12807126783017421E-005,  0.12990189632988428E-005,  0.13013721279939181E-005,  0.13390015381544172E-005,  0.15888056218920581E-005,  0.24631237911666488E-005,],]
input.density_part=[[  0.10000000000000000E+001,  0.10000000000000000E+001,  0.17699995326295348E+004,  0.17699999997521159E+004,  0.17699999999998740E+004,  0.17699999999999998E+004,  0.17700000000000000E+004,  0.17699999999999998E+004,  0.17699999999999998E+004,  0.17699999999999998E+004,  0.17700000000000000E+004,  0.17700000000000603E+004,  0.17700000000082368E+004,  0.17700000008909481E+004,  0.17700000735185765E+004,  0.17700123738902116E+004,  0.17710775639229466E+004,  0.17865578781664963E+004,  0.18501815262078892E+004,  0.18657537553997486E+004,  0.18668399833636963E+004,  0.16404955842344450E+004,  0.16216271544534775E+004,  0.17857908878816424E+004,  0.18966148012085253E+004,  0.19098329190706479E+004,  0.19069652315004410E+004,  0.17888437551174218E+004,  0.17848542772019864E+004,  0.17996108161667582E+004,  0.18411043892605480E+004,  0.18622259231790226E+004,  0.18957388981190154E+004,  0.19157919838257292E+004,  0.19225138402542868E+004,  0.19245813473630178E+004,  0.19246295018066303E+004,  0.19205449180496373E+004,  0.19008063442209900E+004,  0.18708501356194297E+004,  0.18372164661865525E+004,  0.18177143940058368E+004,  0.18112934803390990E+004,  0.16581112457791071E+004,  0.14271919679194161E+004,  0.13753585431828799E+004,  0.12114586632425910E+004,  0.11765793117780461E+004,  0.12857658494196558E+004,  0.12256897824158236E+004,  0.11641459218713849E+004,  0.13335069891877004E+004,  0.14557922101735060E+004,  0.15702629757001782E+004,  0.16676122819811806E+004,  0.17541876013769590E+004,  0.18196174214312673E+004,  0.18559136700288150E+004,  0.18683058345395541E+004,  0.18676689047642951E+004,  0.18548151423011823E+004,  0.18385437046193099E+004,  0.18157190622793842E+004,  0.17864983160694605E+004,  0.16642230845487836E+004,  0.14997591514097624E+004,  0.14200530334338630E+004,  0.14023739166419159E+004,  0.14021181489993742E+004,  0.13695450799291104E+004,  0.12203535981506816E+004,  0.10590373200751299E+004,],]
input.sig_part=[[  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,],]
# Output data is stored here.
output = Object()
output.vlc_grv=[[  0.98603012157266544E-004,  0.68691305877947073E-004,  0.83292769434542974E-001,  0.56552944691864705E-001,  0.37993325979786836E-001,  0.25515157302055284E-001,  0.32033827409032781E-001,  0.27843932209052371E-001,  0.21666312460577802E-001,  0.16506789725824884E-001,  0.12766818029514442E-001,  0.10129256240044984E-001,  0.81453015991418114E-002,  0.67171815539683478E-002,  0.55648130471066008E-002,  0.34607420969547669E-002,  0.14063488593129029E-002,  0.72409657788686773E-003,  0.54028366997876324E-003,  0.48216676451503165E-003,  0.43656470465360046E-003,  0.39182973339325612E-003,  0.36650104084495408E-003,  0.35432114302637579E-003,  0.34810148807543875E-003,  0.33408374404727734E-003,  0.31875471188233626E-003,  0.29061986615985352E-003,  0.31107996883081117E-003,  0.30347970815719613E-003,  0.30697817724473496E-003,  0.30162440416679699E-003,  0.29231253020880468E-003,  0.28121986851494511E-003,  0.27081492080955466E-003,  0.26098123217937341E-003,  0.25105456702234082E-003,  0.24129617006381606E-003,  0.23254090965046313E-003,  0.22343050799646840E-003,  0.21527802704841259E-003,  0.20771798550147760E-003,  0.20206224772245372E-003,  0.20611955617297976E-003,  0.23135259354967424E-003,  0.23772649819866684E-003,  0.28111698715101684E-003,  0.28962411177352765E-003,  0.25396976882671816E-003,  0.30270069393162721E-003,  0.35735339988083139E-003,  0.27242178169860251E-003,  0.25078300807052035E-003,  0.23686360926078313E-003,  0.22609914424589670E-003,  0.21885678395518245E-003,  0.21517323330417514E-003,  0.21110518281121303E-003,  0.20815370114666277E-003,  0.20547510980750277E-003,  0.20368514926526524E-003,  0.20271809137232453E-003,  0.19982999276445083E-003,  0.19911075644119731E-003,  0.20250315399990648E-003,  0.21114411072133951E-003,  0.22050669973799208E-003,  0.22315407845433625E-003,  0.22338634155396833E-003,  0.23003325697531401E-003,  0.28450580897040054E-003,  0.57614340645000851E-003,],]