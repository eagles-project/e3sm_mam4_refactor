# Data at [lat=  -4.98879105856447      , lon=   28.0026147968828      , k=          52 , time step=        1410 ] 
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
input.tair=[[  0.25425190590769432E+003,  0.26867365766471181E+003,  0.27499369167046177E+003,  0.27659787941855774E+003,  0.27105700007007169E+003,  0.26565715589566423E+003,  0.26268081440832935E+003,  0.25773721748988623E+003,  0.24994313884838843E+003,  0.24203196760972332E+003,  0.23621048172803458E+003,  0.23410452519960867E+003,  0.23243743191876754E+003,  0.23004253512817851E+003,  0.22502005024863223E+003,  0.21860515041854060E+003,  0.21164772843483908E+003,  0.20528817043238820E+003,  0.19948780068996274E+003,  0.19465216695449183E+003,  0.18986194473191424E+003,  0.18729749666189190E+003,  0.18703181637614802E+003,  0.18793384559406718E+003,  0.18945194847271159E+003,  0.19179849733389233E+003,  0.19245723263967091E+003,  0.18938678652794172E+003,  0.19088927680573960E+003,  0.19525170540777094E+003,  0.19885478419787083E+003,  0.20302306431689107E+003,  0.20756882509012462E+003,  0.21223398231548208E+003,  0.21685607523984004E+003,  0.22150917504043230E+003,  0.22624788151981227E+003,  0.23088512478386670E+003,  0.23548373907354284E+003,  0.23985832637499681E+003,  0.24413216500862353E+003,  0.24830421141515248E+003,  0.25268168308750148E+003,  0.25687720355592757E+003,  0.26107601949944655E+003,  0.26525661346270061E+003,  0.26763803123463572E+003,  0.26969498264609712E+003,  0.27166610769235700E+003,  0.27364088052166045E+003,  0.27537299883942927E+003,  0.27885096087063289E+003,  0.28157349249790116E+003,  0.28398351917682692E+003,  0.28611848648316055E+003,  0.28791906987207489E+003,  0.28932921413002930E+003,  0.29041733956494483E+003,  0.29131243869083733E+003,  0.29215777048100364E+003,  0.29293508195384510E+003,  0.29364252265878770E+003,  0.29422358191407835E+003,  0.29473703804672158E+003,  0.29502363267774967E+003,  0.29520206861211324E+003,  0.29556264062110802E+003,  0.29607244280795169E+003,  0.29657063692818002E+003,  0.29599099309050143E+003,  0.29488797385232471E+003,  0.29434353631871016E+003,],]
input.pmid=[[  0.12382541305561677E+002,  0.18282923550684863E+002,  0.26994886212093746E+002,  0.39858170362288227E+002,  0.58850914656475304E+002,  0.86893857004050020E+002,  0.12829949082549251E+003,  0.18943524794063876E+003,  0.27970269352932593E+003,  0.41298331550248361E+003,  0.59684493671180769E+003,  0.83774043789417703E+003,  0.11473787230534558E+004,  0.15333938222216857E+004,  0.19996337978167076E+004,  0.25444696509961595E+004,  0.31593251291285251E+004,  0.38366283094427008E+004,  0.45671197939426802E+004,  0.53309561434392535E+004,  0.61015181670369457E+004,  0.68476390229712551E+004,  0.75355335897260984E+004,  0.81946275124734711E+004,  0.88910543149560908E+004,  0.96466673444937041E+004,  0.10466496728675462E+005,  0.11356000142146328E+005,  0.12321099142746600E+005,  0.13368217993210610E+005,  0.14504326553270042E+005,  0.15736987784007928E+005,  0.17074407952158203E+005,  0.18473503838286440E+005,  0.19918075942798609E+005,  0.21463975351308352E+005,  0.23141254216895351E+005,  0.24961078044727059E+005,  0.26935561257428097E+005,  0.29077847535269750E+005,  0.31402197198016987E+005,  0.33924082960874701E+005,  0.36660293473749705E+005,  0.39629043751560472E+005,  0.42850095767885425E+005,  0.46295620874399465E+005,  0.49840531028012527E+005,  0.53368645981549605E+005,  0.56809294208543004E+005,  0.60128432202406751E+005,  0.63311993003465381E+005,  0.66317137730998889E+005,  0.69101791773265897E+005,  0.71625612784052821E+005,  0.73850961206441585E+005,  0.75743833047264721E+005,  0.77302790445671839E+005,  0.78620997685627124E+005,  0.79844100029862268E+005,  0.81038156543039499E+005,  0.82201250345655120E+005,  0.83331494063502178E+005,  0.84427038268033444E+005,  0.85486072930484908E+005,  0.86506827095100627E+005,  0.87487573471257812E+005,  0.88426640136198883E+005,  0.89322421995462981E+005,  0.90173377288449003E+005,  0.90982767070176953E+005,  0.91676597097466714E+005,  0.92113826429928667E+005,],]
input.radius_part=[[  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.15594030812845691E-005,  0.20056132794152794E-005,  0.23161948443659185E-005,  0.25989940148422778E-005,  0.28579105196837105E-005,  0.30739753029533281E-005,  0.32253627919901835E-005,  0.33383526213180592E-005,  0.33581389134597751E-005,  0.26753425037159395E-005,  0.14997211727271614E-005,  0.98363621653527429E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.91219108574494577E-006,  0.91582774041860424E-006,  0.86154063253170507E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.83953472836494377E-006,  0.84765434161746766E-006,  0.90970175040145248E-006,  0.91642596758477417E-006,  0.93274416566217626E-006,  0.93412217708834835E-006,  0.92643151966372540E-006,  0.92307668637052051E-006,  0.92169572878020868E-006,  0.92038568226667546E-006,  0.91880709220856364E-006,  0.91861404316977949E-006,  0.12076993290554310E-005,  0.92598280099315617E-006,  0.92975201883315991E-006,  0.93079478499638410E-006,  0.93192570972915191E-006,  0.10174657818981059E-005,  0.11877601786684681E-005,  0.12204524701335349E-005,  0.14727327050240200E-005,  0.15150111017505913E-005,  0.13362870319434328E-005,  0.15170871010173980E-005,  0.17070063427735205E-005,  0.14025408387417735E-005,  0.12907491455888169E-005,  0.12186476517000325E-005,  0.11629208985423799E-005,  0.11187770465458131E-005,  0.10901063648997656E-005,  0.10698051041079402E-005,  0.10598847119161793E-005,  0.10539528722219727E-005,  0.10543856896166490E-005,  0.10579889362638243E-005,  0.10686170395633708E-005,  0.10753357234857509E-005,  0.11047350092685889E-005,  0.11855637186294895E-005,  0.12701844753773858E-005,  0.13062125627205071E-005,  0.13136327943042542E-005,  0.14313171365864935E-005,  0.16882934057498976E-005,  0.20366094173415824E-005,],]
input.density_part=[[  0.10000000000000000E+001,  0.10000000000000000E+001,  0.17699998116186441E+004,  0.17699999998409724E+004,  0.17699999999998884E+004,  0.17699999999999995E+004,  0.17700000000000000E+004,  0.17700000000000000E+004,  0.17700000000000002E+004,  0.17700000000000000E+004,  0.17700000000000005E+004,  0.17700000000000521E+004,  0.17700000000068624E+004,  0.17700000007347830E+004,  0.17700000608973185E+004,  0.17700111154359631E+004,  0.17710413747026203E+004,  0.17859827633942580E+004,  0.18494452554891554E+004,  0.18649966218811380E+004,  0.18655019752965268E+004,  0.16206506363606559E+004,  0.16051364523957366E+004,  0.17758790365354193E+004,  0.18965171943765147E+004,  0.19075930787743191E+004,  0.19031544738106061E+004,  0.17861603986919172E+004,  0.17858136336342668E+004,  0.18016009366054984E+004,  0.18420386691755427E+004,  0.18725134866185622E+004,  0.19109826045193934E+004,  0.19215292120249537E+004,  0.19239061706148327E+004,  0.19241656144187345E+004,  0.19208605577494357E+004,  0.19099924147579095E+004,  0.13959608043035421E+004,  0.18546592193501951E+004,  0.18268967221356204E+004,  0.18137906272651344E+004,  0.18121453805478732E+004,  0.16442686914145158E+004,  0.14255987641015142E+004,  0.13968352295200880E+004,  0.12112114985239175E+004,  0.11830641168746261E+004,  0.13101011225928214E+004,  0.12452559563966536E+004,  0.11792285337060346E+004,  0.13486975365262019E+004,  0.14718750685566370E+004,  0.15778652504976280E+004,  0.16678483963666326E+004,  0.17523604634224200E+004,  0.18202226145311358E+004,  0.18596831024083385E+004,  0.18737644712314273E+004,  0.18755618306422468E+004,  0.18637024810158944E+004,  0.18481555463671830E+004,  0.18004057244429894E+004,  0.17739505957202859E+004,  0.16999152436130662E+004,  0.15390449029938579E+004,  0.14286380842012120E+004,  0.13938952827365345E+004,  0.13898167463089167E+004,  0.12996041213343447E+004,  0.11725482958458274E+004,  0.10949833907887360E+004,],]
input.sig_part=[[  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,  0.18000000000000000E+001,],]
# Output data is stored here.
output = Object()
output.vlc_grv=[[  0.27805427065870136E-003,  0.19366379870475466E-003,  0.23502395129033318E+000,  0.15979666721307276E+000,  0.10730874202498938E+000,  0.72125666124597310E-001,  0.91398959732902116E-001,  0.80413522851213615E-001,  0.63939635017874127E-001,  0.50493365838923547E-001,  0.41287164767636740E-001,  0.35242105099026880E-001,  0.30926231837528954E-001,  0.28090506816952646E-001,  0.25577290903276276E-001,  0.16325413865521460E-001,  0.59931027278121137E-002,  0.29172288021047537E-002,  0.21950534441750550E-002,  0.20740884028896374E-002,  0.19805207157246323E-002,  0.19062063467196783E-002,  0.18470673858068389E-002,  0.17923163512776814E-002,  0.17818796931417644E-002,  0.17420048632192247E-002,  0.16988705362237615E-002,  0.16049697969851401E-002,  0.17807893046018583E-002,  0.17670495592102950E-002,  0.18166644743835341E-002,  0.18015392153965957E-002,  0.17611157200737720E-002,  0.17120820525121467E-002,  0.16660958035995999E-002,  0.16208149855962927E-002,  0.15734744777708230E-002,  0.15270641348682668E-002,  0.18259288801847793E-002,  0.14386787017250208E-002,  0.13981342801934558E-002,  0.13630606659461512E-002,  0.13374531865043380E-002,  0.14067325475277768E-002,  0.16109072364041482E-002,  0.16349039006114860E-002,  0.20169942502867629E-002,  0.20622565332819700E-002,  0.17731886882500006E-002,  0.21404348838093118E-002,  0.25352112681814526E-002,  0.19520773434806168E-002,  0.17950958374400949E-002,  0.17065838217049988E-002,  0.16351479734029327E-002,  0.15839628018566273E-002,  0.15570182724488426E-002,  0.15280862283061827E-002,  0.15073776655343903E-002,  0.14881170462624879E-002,  0.14759149858290597E-002,  0.14697155064501036E-002,  0.14567956213751914E-002,  0.14502550162930996E-002,  0.14628487138493461E-002,  0.15186779680441828E-002,  0.16110170030545457E-002,  0.16575211440939782E-002,  0.16683997794751616E-002,  0.18480116764710837E-002,  0.23125913951274140E-002,  0.31298316312454218E-002,],]