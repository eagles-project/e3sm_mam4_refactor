# Data at [lat=  -4.98879105856447      , lon=   28.0026147968828      , k=          52 , time step=        1389 ] 
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
input.tair=[[  0.25255095988639337E+003,  0.26107404234460671E+003,  0.26348899752764191E+003,  0.26106884668874710E+003,  0.26307855343404026E+003,  0.26540770958074438E+003,  0.26391778285816054E+003,  0.25882422327161731E+003,  0.24870197818809211E+003,  0.24018422563156011E+003,  0.23565178302108231E+003,  0.23551749575176052E+003,  0.23428151510000140E+003,  0.23131651710640051E+003,  0.22422576665846250E+003,  0.21608550184521152E+003,  0.21021292266057139E+003,  0.20629645031375955E+003,  0.20172601341021070E+003,  0.19722344949792395E+003,  0.19231887614256348E+003,  0.18924055680923709E+003,  0.18822281827066104E+003,  0.18828970532424825E+003,  0.18855128078255058E+003,  0.18878257596502124E+003,  0.18667113395685715E+003,  0.18799907344425731E+003,  0.19134279695252144E+003,  0.19513748524335679E+003,  0.19934710222395117E+003,  0.20353795447376979E+003,  0.20790319830536464E+003,  0.21223713737453176E+003,  0.21667152718647489E+003,  0.22121732000831378E+003,  0.22588248922357610E+003,  0.23048719793169269E+003,  0.23515835867395833E+003,  0.23977875002299015E+003,  0.24422700284787990E+003,  0.24822946396596251E+003,  0.25246309326692403E+003,  0.25701704885815030E+003,  0.26179574992192454E+003,  0.26551925040793628E+003,  0.26790748020975809E+003,  0.27032406351803269E+003,  0.27139982300974339E+003,  0.27319351071038614E+003,  0.27526240263210735E+003,  0.27899983325741732E+003,  0.28167977742252850E+003,  0.28401091565803796E+003,  0.28602366791747556E+003,  0.28760566081331478E+003,  0.28882984279587345E+003,  0.28979920941997568E+003,  0.29064136959715904E+003,  0.29141697872997099E+003,  0.29219013720989432E+003,  0.29289607515061869E+003,  0.29348934445837426E+003,  0.29398454640435665E+003,  0.29443691780866584E+003,  0.29486179442743378E+003,  0.29529757451005912E+003,  0.29573751046014809E+003,  0.29564280560384321E+003,  0.29419658512940970E+003,  0.29338902498761416E+003,  0.29326896336311432E+003,],]
input.pmid=[[  0.12382541305561677E+002,  0.18282923550684863E+002,  0.26994886212093746E+002,  0.39858170362288227E+002,  0.58850914656475304E+002,  0.86893857004050020E+002,  0.12829949082549251E+003,  0.18943524794063876E+003,  0.27970269352932593E+003,  0.41298331550248361E+003,  0.59684493671180769E+003,  0.83774043789417703E+003,  0.11473787230534558E+004,  0.15333938222216857E+004,  0.19996337978167076E+004,  0.25444696509961595E+004,  0.31593251291285251E+004,  0.38366283094427008E+004,  0.45671197939426802E+004,  0.53309561434392535E+004,  0.61015181670369457E+004,  0.68476390229712551E+004,  0.75355335897260984E+004,  0.81946275124734711E+004,  0.88910543149560908E+004,  0.96466673444937041E+004,  0.10466496728675462E+005,  0.11356000142146328E+005,  0.12321099142746600E+005,  0.13368217993210610E+005,  0.14504326553270042E+005,  0.15736987784007928E+005,  0.17074407952158203E+005,  0.18474053418759373E+005,  0.19919998054664902E+005,  0.21467613304312126E+005,  0.23146753833705159E+005,  0.24968597540405488E+005,  0.26945272293188202E+005,  0.29089936360655087E+005,  0.31416865889929850E+005,  0.33941550770955429E+005,  0.36680798288504666E+005,  0.39652843673818730E+005,  0.42877470834413325E+005,  0.46326820235122439E+005,  0.49875664994468709E+005,  0.53407695912173796E+005,  0.56852163020775850E+005,  0.60174985028292002E+005,  0.63362079361088741E+005,  0.66370559590747958E+005,  0.69158304406262570E+005,  0.71684926684208127E+005,  0.73912745089781689E+005,  0.75807717887809456E+005,  0.77368405621417973E+005,  0.78688075980189780E+005,  0.79912535882699216E+005,  0.81107917715411153E+005,  0.82272302471522067E+005,  0.83403800681631241E+005,  0.84500560865017280E+005,  0.85560770983195776E+005,  0.86582658113992933E+005,  0.87564493050482721E+005,  0.88504602014964606E+005,  0.89401378130572251E+005,  0.90253277924664828E+005,  0.91063510657475214E+005,  0.91757989066597074E+005,  0.92195606578348583E+005,],]
input.radius_part=[[  0.21709496539507993E-007,  0.22566096639251814E-007,  0.26028503904743299E-007,  0.29239690552391464E-007,  0.31167849136310405E-007,  0.35604855258690754E-007,  0.40155826071256915E-007,  0.42589798968869089E-007,  0.44649285614410032E-007,  0.46419358557947579E-007,  0.48158902218786721E-007,  0.49788397373341782E-007,  0.51069000094988676E-007,  0.52179339364855610E-007,  0.54760258776648729E-007,  0.59533812693087371E-007,  0.64282046486810310E-007,  0.66844684512465580E-007,  0.67974416829490588E-007,  0.69189206732121961E-007,  0.69500095115741527E-007,  0.69642718832189588E-007,  0.69642719151116048E-007,  0.69642718880887047E-007,  0.69642718832189588E-007,  0.69599459008020050E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642718832189588E-007,  0.69642719640992921E-007,  0.69642719754538852E-007,  0.69642720286869665E-007,  0.69642720134984561E-007,  0.69642720697098905E-007,  0.69642722785073276E-007,  0.69642722699223535E-007,  0.69642720025672616E-007,  0.69642719500210261E-007,  0.69642719266998639E-007,  0.69642719156091917E-007,  0.69642719145331554E-007,  0.69642719161735168E-007,  0.69642719169207772E-007,  0.69642719178027026E-007,  0.69642719182511708E-007,  0.69642719171618997E-007,  0.69642719214498120E-007,  0.69642719244191231E-007,  0.69642719278189082E-007,  0.69642719382321619E-007,  0.69642719610096296E-007,  0.69642719839893753E-007,  0.69642719927515487E-007,  0.69642719919910799E-007,  0.69642720360768869E-007,  0.69642722020778176E-007,  0.69642722516418221E-007,],]
input.density_part=[[  0.10465307489579268E+004,  0.10463926832310920E+004,  0.10464948281552740E+004,  0.10462222653211136E+004,  0.10451596270269406E+004,  0.10451622827343479E+004,  0.10453141138351671E+004,  0.10449932234776454E+004,  0.10445832756175421E+004,  0.10446923045820165E+004,  0.10436794787120061E+004,  0.10421509369690896E+004,  0.10406345773277612E+004,  0.10370761115556204E+004,  0.10334114855027181E+004,  0.10378502872163551E+004,  0.10501243428060368E+004,  0.10662786984004038E+004,  0.10682478329190660E+004,  0.10628412134242121E+004,  0.10575670228385977E+004,  0.10565185730821374E+004,  0.10565213992438698E+004,  0.10581390353953145E+004,  0.10591030750696577E+004,  0.10592613098516131E+004,  0.10635103056386358E+004,  0.10645524631837227E+004,  0.10638427191385617E+004,  0.10633647669136124E+004,  0.10620991789892792E+004,  0.10618961855847331E+004,  0.10615510535245448E+004,  0.10603148019609889E+004,  0.10599470718977689E+004,  0.10600795650703419E+004,  0.10595317501288011E+004,  0.10582228421259204E+004,  0.10577940526724931E+004,  0.10574114912975360E+004,  0.10567485786105360E+004,  0.10575970592023164E+004,  0.10575164882573595E+004,  0.10573278022336963E+004,  0.10571265214184421E+004,  0.10577432925233452E+004,  0.10563936542389342E+004,  0.10495329467633951E+004,  0.10473638633511773E+004,  0.10448041548980245E+004,  0.10437154694681524E+004,  0.10432856533391166E+004,  0.10432711863754716E+004,  0.10432811707449218E+004,  0.10436093557385886E+004,  0.10434712340786898E+004,  0.10438399458542572E+004,  0.10442847592555431E+004,  0.10444338835690476E+004,  0.10445936191195303E+004,  0.10451102232632109E+004,  0.10452886333264678E+004,  0.10454782456087107E+004,  0.10472110284579931E+004,  0.10498724156791011E+004,  0.10528798649419844E+004,  0.10565114766303932E+004,  0.10589526200593730E+004,  0.10595003187940006E+004,  0.10615889065343406E+004,  0.10638420762853777E+004,  0.10640700485397242E+004,],]
input.sig_part=[[  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,  0.16000000238418579E+001,],]
# Output data is stored here.
output = Object()
output.vlc_grv=[[  0.24967905281314869E-002,  0.17869212646175535E-002,  0.14025208410688957E-002,  0.10619129390439643E-002,  0.76881708125028220E-003,  0.59748518947602680E-003,  0.45520734894638476E-003,  0.32375933431032483E-003,  0.22529814990633959E-003,  0.15596758694967492E-003,  0.11085066162864405E-003,  0.81556391826228652E-004,  0.60881245969912246E-004,  0.46143901114896441E-004,  0.36495141876673904E-004,  0.30821793093523142E-004,  0.26837721494613532E-004,  0.23196814181683379E-004,  0.19706362497426171E-004,  0.16979264981953326E-004,  0.14708451119208602E-004,  0.13069824283100178E-004,  0.11887079589728439E-004,  0.10984461530346301E-004,  0.10173811237889050E-004,  0.94117097720483574E-005,  0.87077514919696539E-005,  0.80932785007047319E-005,  0.75451524145367476E-005,  0.70429770451035682E-005,  0.65752502333520266E-005,  0.61445970021227117E-005,  0.57435592597749810E-005,  0.53779032968378832E-005,  0.50567305694575092E-005,  0.47604207791930467E-005,  0.44777133324330117E-005,  0.42068812214600929E-005,  0.39550048272506634E-005,  0.37172293609689904E-005,  0.34913782534385157E-005,  0.32816350173920301E-005,  0.30827998596091259E-005,  0.28970542867267282E-005,  0.27234203711170773E-005,  0.25616818224194269E-005,  0.24099361742457062E-005,  0.22668988363511158E-005,  0.21506079603380355E-005,  0.20515448954070812E-005,  0.19693209605537916E-005,  0.19028111355851247E-005,  0.18455861456114773E-005,  0.17973386562478075E-005,  0.17579231330501684E-005,  0.17253056062156152E-005,  0.17002901354553545E-005,  0.16800284488171359E-005,  0.16613296040794734E-005,  0.16435770008996369E-005,  0.16273529067299138E-005,  0.16114678976914824E-005,  0.15963993467341171E-005,  0.15844603922015074E-005,  0.15747041265952768E-005,  0.15662258480545611E-005,  0.15594445160617499E-005,  0.15516581805546289E-005,  0.15413704535145961E-005,  0.15328264707413957E-005,  0.15267127513574413E-005,  0.15215559026335947E-005,],]