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
input.moment=[[         3],]
input.ncol=[[         4],]
input.pcols=[[         4],]
input.nver=[[        72],]
input.radius_max=[[  0.50000000000000002E-004],]
input.tair=[[  0.25448531923775388E+003,  0.26917507711550218E+003,  0.27527995484790199E+003,  0.27646320368299450E+003,  0.27172252748854567E+003,  0.26670929531280802E+003,  0.26359229089452674E+003,  0.25827808388259535E+003,  0.24999793370358381E+003,  0.24194112389804837E+003,  0.23612810366470794E+003,  0.23404715009753059E+003,  0.23239424896390403E+003,  0.23000030395266378E+003,  0.22488400259528524E+003,  0.21831966249088592E+003,  0.21125272062882010E+003,  0.20500071905936795E+003,  0.19934742129278064E+003,  0.19459643082999108E+003,  0.18991880780124160E+003,  0.18745431489477437E+003,  0.18720657966674224E+003,  0.18808602612324276E+003,  0.18955005485861017E+003,  0.19189372140924138E+003,  0.19246292185821773E+003,  0.18942419619735230E+003,  0.19092052060648689E+003,  0.19530521835040020E+003,  0.19895372900264405E+003,  0.20312470143001036E+003,  0.20763762754290809E+003,  0.21223112923394919E+003,  0.21679596350125865E+003,  0.22145864194805193E+003,  0.22621343845600012E+003,  0.23091885265942648E+003,  0.23553714421331264E+003,  0.24008905161479922E+003,  0.24440276511803674E+003,  0.24854858448827110E+003,  0.25276771067008607E+003,  0.25693665711977485E+003,  0.26103863891854525E+003,  0.26516762471027391E+003,  0.26765853652637441E+003,  0.26962471703376048E+003,  0.27146797654160008E+003,  0.27352337998539548E+003,  0.27532258695030811E+003,  0.27878137396870255E+003,  0.28147124890647115E+003,  0.28392761488813881E+003,  0.28610158929374558E+003,  0.28792231667936193E+003,  0.28932691344681268E+003,  0.29039396193302127E+003,  0.29128334984127679E+003,  0.29213778156204290E+003,  0.29291487406086998E+003,  0.29361045287862811E+003,  0.29419459517035313E+003,  0.29470365594657596E+003,  0.29494991964494443E+003,  0.29517145660896813E+003,  0.29563598263525000E+003,  0.29621059631874215E+003,  0.29674700642431537E+003,  0.29684167922222639E+003,  0.29552990453632617E+003,  0.29443385975310053E+003,],]
input.pmid=[[  0.12382541305561677E+002,  0.18282923550684863E+002,  0.26994886212093746E+002,  0.39858170362288227E+002,  0.58850914656475304E+002,  0.86893857004050020E+002,  0.12829949082549251E+003,  0.18943524794063876E+003,  0.27970269352932593E+003,  0.41298331550248361E+003,  0.59684493671180769E+003,  0.83774043789417703E+003,  0.11473787230534558E+004,  0.15333938222216857E+004,  0.19996337978167076E+004,  0.25444696509961595E+004,  0.31593251291285251E+004,  0.38366283094427008E+004,  0.45671197939426802E+004,  0.53309561434392535E+004,  0.61015181670369457E+004,  0.68476390229712551E+004,  0.75355335897260984E+004,  0.81946275124734711E+004,  0.88910543149560908E+004,  0.96466673444937041E+004,  0.10466496728675462E+005,  0.11356000142146328E+005,  0.12321099142746600E+005,  0.13368217993210610E+005,  0.14504326553270042E+005,  0.15736987784007928E+005,  0.17074407952158203E+005,  0.18473058038959040E+005,  0.19916516796832941E+005,  0.21461024378553546E+005,  0.23136793131402825E+005,  0.24954978508297114E+005,  0.26927684024841717E+005,  0.29068041527794463E+005,  0.31390298498353026E+005,  0.33909913720210810E+005,  0.36643660727361421E+005,  0.39609738135991014E+005,  0.42827890127680330E+005,  0.46270313108365823E+005,  0.49812031655515690E+005,  0.53336970123936953E+005,  0.56774520615591464E+005,  0.60090670274007964E+005,  0.63271364804975536E+005,  0.66273803896706071E+005,  0.69055950818173136E+005,  0.71577499544778417E+005,  0.73800844408276869E+005,  0.75692012030418933E+005,  0.77249565844989003E+005,  0.78566586257143019E+005,  0.79788587400379591E+005,  0.80981568863429944E+005,  0.82143615492377197E+005,  0.83272841612557997E+005,  0.84367399460072804E+005,  0.85425480636368346E+005,  0.86445315780774035E+005,  0.87425179157101345E+005,  0.88363400347253671E+005,  0.89258375702591526E+005,  0.90108564851256990E+005,  0.90917270862283287E+005,  0.91610574946241948E+005,  0.92047489402039384E+005,],]
input.radius_part=[[  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,  0.50000000000000004E-005,],]
input.density_part=[[  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,  0.10000000000000000E+004,],]
input.sig_part=[[  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,  0.14600000000000000E+001,],]
# Output data is stored here.
output = Object()
output.vlc_grv=[[  0.81917704603348040E+000,  0.57160357696588882E+000,  0.39261014488847507E+000,  0.26765419748603403E+000,  0.18099736109130563E+000,  0.12278199053442346E+000,  0.84029484288705214E-001,  0.57804015151562263E-001,  0.40158627709669215E-001,  0.28565789062082823E-001,  0.21399730141163094E-001,  0.16997557826300354E-001,  0.14186547701907888E-001,  0.12384696779562950E-001,  0.11270011574066638E-001,  0.10608626740085768E-001,  0.10239091363106063E-001,  0.10029920149960448E-001,  0.99251081084367684E-002,  0.98778547568686495E-002,  0.98934205723245506E-002,  0.98718339952306149E-002,  0.97941323438681602E-002,  0.96887607549337451E-002,  0.95676958058184582E-002,  0.94181123585127605E-002,  0.93422201137748472E-002,  0.94113437690260815E-002,  0.93064823382264950E-002,  0.90958782474332561E-002,  0.89222743771781403E-002,  0.87382818089809747E-002,  0.85516296040470761E-002,  0.83725544058592654E-002,  0.82041112093598750E-002,  0.80407632517698509E-002,  0.78821779355050386E-002,  0.77319908478612580E-002,  0.75905813691457867E-002,  0.74568018454146861E-002,  0.73341788823325310E-002,  0.72202144651034956E-002,  0.71089578619671469E-002,  0.70029798094299036E-002,  0.69023089615386007E-002,  0.68049194142301116E-002,  0.67435408814917556E-002,  0.66947910110593863E-002,  0.66504065378378335E-002,  0.66037634728686127E-002,  0.65634987406763120E-002,  0.64943943191452065E-002,  0.64413127876247487E-002,  0.63939917357958582E-002,  0.63529721501824170E-002,  0.63191647696941447E-002,  0.62933064735155213E-002,  0.62736941685835913E-002,  0.62573303884713724E-002,  0.62417144762373375E-002,  0.62275264650412500E-002,  0.62148109512387816E-002,  0.62040316283782137E-002,  0.61945764734478550E-002,  0.61894591575337348E-002,  0.61848160110131133E-002,  0.61763370572304492E-002,  0.61661837211006822E-002,  0.61567388664841762E-002,  0.61544128825976518E-002,  0.61747318320218673E-002,  0.61920057926661441E-002,],]