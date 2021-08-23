import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

np.seterr(divide='ignore')

a = 20
b = 78
L = 105

#zeros = [[0.05302558037807864, 0.10766084030369445, 0.16202764889641477, 0.216302979499391, 0.2705371060280038, 0.3247493649044931, 0.37894869199949327, 0.4331397608818642, 0.487325244576556, 0.5415067790951656, 0.5956854204067052, 0.6498618798980598, 0.7040366539963723, 0.7582100994845619, 0.8123824792850834, 0.8665539913646428, 0.9207247875976772, 0.9748949864210473, 1.0290646815465996, 1.0832339480831659, 1.1374028469315154, 1.1915714279872687, 1.2457397325175266, 1.299907794947415, 1.3540756442201864, 1.4082433048456615, 1.4624107977175065, 1.5165781407558478, 1.5707453494178198, 1.6249124371061394], [0.05740913834160282, 0.11031590952507575, 0.1638955208133877, 0.2177351020219049, 0.2716955132827949, 0.3257207720907756, 0.37978456739633454, 0.43387303245350595, 0.48797820486255605, 0.5420952003608704, 0.5962208620925039, 0.650353061473256, 0.7044903122364792, 0.758631545751189, 0.8127759742825448, 0.8669230045707884, 0.921072181410012, 0.9752231497626213, 1.029375628683828, 1.0835293929763405, 1.1376842600335968, 1.191840080236322, 1.2459967298333898, 1.300154105587299, 1.3543121206974666, 1.4084707016560896, 1.462629785799319, 1.5167893193805526, 1.5709492560392615, 1.6251095555755262], [0.06841073580473017, 0.11790099943697449, 0.1694000929194389, 0.22199627720783213, 0.27515516705498666, 0.3286269569353504, 0.38228755663593644, 0.43606994238121044, 0.48993514878548894, 0.543859108341318, 0.5978262007762416, 0.6518258658147359, 0.705850716712304, 0.7598954357771306, 0.8139560995829912, 0.868029751337456, 0.9221141211507546, 0.9762074379231439, 1.0303082997204238, 1.0844155825097994, 1.1385283746516572, 1.1926459290628304, 1.2467676277357336, 1.3008929550460668, 1.3550214774173657, 1.4091528276390592, 1.463286692649097, 1.5174228039147546, 1.5715609297938553, 1.62570086941442], [0.08247156972738798, 0.12938618237602922, 0.17824504002807104, 0.22898109555266666, 0.28087000181309557, 0.3334443407780454, 0.386444108353404, 0.4397220086384786, 0.4931904158163463, 0.5467945343848483, 0.6004985436876339, 0.6542781173325819, 0.7081161886564592, 0.7620004461238923, 0.8159217927914741, 0.8698733654811863, 0.9238498910499622, 0.9778472524289602, 1.0318621890396256, 1.0858920855599652, 1.1399348201354556, 1.1939886534330835, 1.2480521462910277, 1.3021240977331103, 1.3562034977127773, 1.4102894906646843, 1.4643813470903777, 1.5184784411857455, 1.5725802330653054, 1.6266862545190786], [0.09744025930592673, 0.14340052866857, 0.1899310506113393, 0.23850882659578196, 0.2887635325459192, 0.3401348363194361, 0.39223265596909734, 0.44481580776769886, 0.4977350503965873, 0.5508951848714021, 0.6042332873350336, 0.6577063405569306, 0.7112840361328108, 0.7649444471278634, 0.8186713388025004, 0.8724524446033008, 0.9262783292384132, 0.9801416195764336, 1.034036472331197, 1.0879581979852357, 1.1419029901321247, 1.1958677274031104, 1.2498498262967894, 1.303847130314279, 1.3578578253803668, 1.4118803745712034, 1.4659134672026932, 1.5199559787293968, 1.5740069388711733, 1.6280655060640097], [0.11248638134122295, 0.1586715543178695, 0.20376233422056608, 0.2503025932427223, 0.2987230560396181, 0.3486455538578608, 0.39962445166568716, 0.45133380653785077, 0.503557433736567, 0.5561528896882842, 0.6090244313175168, 0.6621059807525411, 0.7153507110134838, 0.7687246162475573, 0.8222024529037087, 0.8757651119369099, 0.9293978744339053, 0.9830892261453567, 1.0368300342514083, 1.0906129642951221, 1.1444320597444548, 1.198282433839355, 1.252160040351506, 1.306061500721041, 1.3599839720689482, 1.4139250452697991, 1.4678826653990102, 1.5218550690418668, 1.5758407344433807, 1.6298383415334803], [0.1273919665478332, 0.17435500508573568, 0.21894318892427386, 0.2639595975138969, 0.31058139933717, 0.3589038191004235, 0.40858309437238255, 0.4592549517849721, 0.510643961008548, 0.5625581522670602, 0.6148649874644915, 0.6674717008669012, 0.7203120240869428, 0.7733375950411581, 0.8265123966141878, 0.8798091029029662, 0.9332066314440598, 0.986688469584898, 1.0402415073298636, 1.093855208147927, 1.147521009683537, 1.2012318837328473, 1.2549820084409873, 1.308766520831764, 1.3625813276793373, 1.4164229593352629, 1.4702884555703135, 1.524175275557955, 1.5780812262614663, 1.632004404983989], [0.1421340106392586, 0.19004935615642518, 0.23477339395896973, 0.2789586832858553, 0.3240872861448627, 0.37080321374348324, 0.41906058551767816, 0.46855418237296664, 0.5189794563243123, 0.5701006868062326, 0.6217474407976348, 0.6737977379845761, 0.7261634118465511, 0.77877968741993, 0.8315981255378513, 0.8845818762891131, 0.9377024556246906, 0.9909375228737823, 1.0442693223183859, 1.0976835718480047, 1.1511686575652722, 1.2047150412300356, 1.2583148181340575, 1.3119613829033265, 1.365649173835035, 1.4193734751306921, 1.4731302623352864, 1.526916080391146, 1.5807279465680377, 1.6345632725580737], [0.1567321336831411, 0.20562006822466378, 0.25079439602090986, 0.29475145249668633, 0.3388877454958393, 0.38417952294352103, 0.43098694166646717, 0.47919939599890415, 0.5285467431294443, 0.5787697221463958, 0.6296641786879656, 0.6810782872940595, 0.7329002427263638, 0.7850470946862402, 0.8374564677813215, 0.8900807497796659, 0.9428830564508422, 0.9958344143340989, 1.048911770275583, 1.1020965651004675, 1.1553736965568893, 1.2087307545470911, 1.2621574495249877, 1.315645179858834, 1.3691867005125138, 1.422775866535076, 1.4764074324501795, 1.530076893877137, 1.5837803613979706, 1.6375144592781876], [0.1712090105136812, 0.22104291763301306, 0.26678311246205366, 0.31089288903671414, 0.3545620313150986, 0.39878782011946323, 0.4442515969210815, 0.49114373217384294, 0.5393243352293513, 0.5885536534294308, 0.6386077261593192, 0.6893078512653247, 0.740518138943503, 0.7921361784353904, 0.8440843298588003, 0.8963030590302328, 0.9487461202370847, 1.0013771227101613, 1.054167076744661, 1.1070926233568978, 1.1601347416637855, 1.2132777930013054, 1.2665088051080313, 1.3198169295797735, 1.3731930259418037, 1.4266293393664071, 1.4801192484498225, 1.5336570659767037, 1.5872378801675324, 1.6408574271564071], [0.1855833469828302, 0.23632679114403435, 0.2826603024260202, 0.32710792137316197, 0.37071057880513236, 0.41430780194361205, 0.45868160725123786, 0.5043113653890567, 0.5512806700222785, 0.599438289481553, 0.6485704781320762, 0.6984814323239379, 0.7490132666634699, 0.8000437332607608, 0.8514789235136896, 0.903246338213089, 0.9552894516488805, 1.0075636877664114, 1.0600334884501836, 1.1126701761755922, 1.1654503840209396, 1.2183548904651054, 1.2713677448280454, 1.3244756033953735, 1.377667219935678, 1.4309330506162754, 1.4842649446166807, 1.5376558996163, 1.5910998668864338, 1.644591594669166], [0.19986984544095143, 0.2514868208242336, 0.2984099754616655, 0.34326626183784814, 0.38703808408978213, 0.4303960811790128, 0.4740337960936865, 0.5185786943997984, 0.5643632591342397, 0.6114025100880559, 0.6595433695298591, 0.7085943361036009, 0.7583825014380824, 0.8087672328227222, 0.8596379995894758, 0.9109085171329391, 0.9625111326299756, 1.014392336344622, 1.0665093728956012, 1.1188277260798085, 1.1713192536448354, 1.2239607956253196, 1.2767331266193476, 1.3296201590184233, 1.382608330839093, 1.4356861306575504, 1.4888437253989888, 1.5420726660606243, 1.5953656530648301, 1.648716347656927], [0.21408012748075317, 0.266537273509463, 0.31403735203279054, 0.35932184677512663, 0.4033741270250389, 0.44676021807171995, 0.49002096344216955, 0.5337612430421577, 0.5784831030841868, 0.6244099479784319, 0.671512579874859, 0.7196411622174487, 0.7686232906014062, 0.8183049771126912, 0.8685600604905884, 0.9192881241558579, 0.9704096950748801, 1.0218616229986037, 1.0735933310098258, 1.1255639383785891, 1.1777400911604305, 1.230094329514105, 1.2826038528250248, 1.3352495782452178, 1.38801541635449, 1.4408877086161223, 1.4938547864327787, 1.5469066224406347, 1.600034552415683, 1.6532310517262716], [0.22822352763783668, 0.2814903108744866, 0.32955268312410935, 0.375265796528306, 0.4196414114474255, 0.46320309439165497, 0.5063683738370786, 0.549623096985803, 0.5935003098954721, 0.6383963197573608, 0.6844531217087572, 0.7316132864201553, 0.7797328878708588, 0.828656002592081, 0.8782444946170633, 0.9283844724861121, 0.9789842988043167, 1.0299705824372052, 1.0812843222997859, 1.132877742226997, 1.1847118290127923, 1.2367544508380288, 1.288978922957836, 1.3413629097993331, 1.3938875785488096, 1.446536941157804, 1.499297338365089, 1.5521570315905857, 1.605105877472317, 1.658135066250698], [0.2423076660814645, 0.2963561330769698, 0.3449662858307125, 0.3911022462796289, 0.4358136430739078, 0.4796184515317208, 0.5228653576308095, 0.5659154993899279, 0.6092235338385896, 0.6532556351857703, 0.6983186346028812, 0.7444938906009897, 0.7917064184429536, 0.8398195029420376, 0.8886915242948462, 0.938197784944985, 0.9882348969785425, 1.0387188871187147, 1.0895818004279145, 1.140768442632409, 1.1922336825187547, 1.2439403296337563, 1.295857493328607, 1.3479593177811684, 1.400224003419748, 1.4526330449906697, 1.5051706337143769, 1.5578231843834902, 1.6105789582709948, 1.6634277600871739], [0.2563388542358683, 0.31114333375246317, 0.360287378240645, 0.40683870860503957, 0.45188670080928733, 0.49596048296645767, 0.5393814943937814, 0.5824255143593348, 0.6254335535072951, 0.6688341606883826, 0.7130289765002551, 0.7582497729417295, 0.8045330056283974, 0.8517933352333545, 0.8999017692674741, 0.9487291714008808, 0.9981623531386576, 1.048106995678179, 1.0984858539309306, 1.149235841845329, 1.20030525065211, 1.2516514296566525, 1.3032389440710375, 1.3550381362269321, 1.407024005462526, 1.4591753334439153, 1.5114739970626474, 1.5639044247952258, 1.6164531632817585, 1.669108529155759], [0.2703223860785714, 0.3258592180505185, 0.3755239839789204, 0.42248290050408255, 0.467864328938827, 0.512214657953326, 0.5558508819862268, 0.5990048402862147, 0.6419239255212267, 0.6849426397455193, 0.7284615714901469, 0.7728205651512298, 0.818189206262796, 0.8645709863975455, 0.9118750894501508, 0.9599803021737952, 1.0087684413802624, 1.0581362634057845, 1.1079973405452082, 1.1582803659262337, 1.2089266253909723, 1.2598875995232017, 1.3111229540017164, 1.362598930308977, 1.4142870777273806, 1.4661632575478238, 1.5182068589258426, 1.5704001779800785, 1.622727922810718, 1.6751768160580562], [0.2842627515342311, 0.3405100530561104, 0.39068307637276656, 0.4380418994103735, 0.48375237779439995, 0.5283797527601156, 0.572246709068884, 0.6155684432893886, 0.658535246320798, 0.7013872871749398, 0.7444556518212941, 0.788109190405095, 0.8326298319595633, 0.8781363092678748, 0.9246081988788903, 0.9719525096037548, 1.0200556047168492, 1.068808959757701, 1.1181179918721107, 1.1679031871556391, 1.2180985062736855, 1.268649171746064, 1.3195095834234842, 1.3706415646173908, 1.4220129478814807, 1.473596452087462, 1.5253687937027953, 1.5773099816859404, 1.6294027551411874, 1.6816321319475875], [0.29816379597243725, 0.35510126091795924, 0.40577074557004694, 0.4535220164439253, 0.49955683647983373, 0.5444592535402389, 0.5885614056604235, 0.6320759314299087, 0.6751656587704649, 0.7180040975045077, 0.7608339380325969, 0.80398085982047, 0.847778505804165, 0.8924558079358091, 0.9380904395673736, 0.9846449044095094, 1.0320262570143055, 1.080128092583394, 1.128850443928589, 1.1781063233918265, 1.2278223125895367, 1.2779370669870853, 1.3283993642558873, 1.37916627870777, 1.4302016407216507, 1.4814747871182716, 1.5329595621569294, 1.5846335213918397, 1.636477295724805, 1.6884740809045713], [0.31202884133359715, 0.3696375686183643, 0.42079234311507924, 0.46892885606318996, 0.5152832741361234, 0.5604579003888607, 0.6047954003325948, 0.6485119827902606, 0.6917606579747695, 0.7346776897095754, 0.777433518341139, 0.8202767581310674, 0.8635234594192446, 0.9074698181865744, 0.9522973490043793, 0.9980509772856402, 1.0446812915559858, 1.0920968661040296, 1.140198112540213, 1.1888926782526628, 1.2381002777836252, 1.2877528973615098, 1.337793395361815, 1.3881737694840508, 1.4388535473664192, 1.4897984253909822, 1.5409791588963277, 1.592370669593157, 1.6439513297875117, 1.6957023871136208], [0.3258607796777946, 0.3841231255390591, 0.43575259982171455, 0.48426740830621345, 0.5309367381881868, 0.5763804835274006, 0.6209519164318946, 0.6648731001961782, 0.7082959114182894, 0.751339760674584, 0.7941293620817563, 0.8368396743522534, 0.8797248725146138, 0.923086555417165, 0.9671827067469907, 1.012153262304709, 1.0580173434804117, 1.1047174973700817, 1.1521647734871587, 1.2002659563252847, 1.2489354965995594, 1.298099056169709, 1.3476934380385284, 1.3976652787197834, 1.4479695008466353, 1.4985678859447322, 1.549427865289567, 1.600521529686086, 1.6518248287505308, 1.7033169251834765], [0.3396621463666888, 0.3985615967559572, 0.45065572109313135, 0.4995421352371885, 0.5465217813538389, 0.5922314888815555, 0.6370348924102918, 0.6811605860652752, 0.7247630179277735, 0.7679560168782407, 0.8108404963857507, 0.8535389733880112, 0.8962346664576597, 0.9391847213751943, 0.9826718143820062, 1.0269163257791207, 1.0720223760990146, 1.1179890072689165, 1.1647536218975922, 1.2122303393633336, 1.2603318712667673, 1.308978769900735, 1.3581020013092933, 1.4076426816157688, 1.4575508567523585, 1.5077841137054386, 1.5583063080993134, 1.6090864848674695, 1.6600979908935725, 1.7113177539923248], [0.3534351781019158, 0.41295623782300445, 0.465505464414199, 0.5147570446304037, 0.5620425161305531, 0.6080150284439494, 0.6530482503412032, 0.6973773819216518, 0.74116081966907, 0.7845121534974625, 0.8275221488605079, 0.8702828297982057, 0.9129204767067677, 0.9556276789640862, 0.9986607732924381, 1.042279911707972, 1.0866695948046656, 1.1319035768844148, 1.1779654870144298, 1.2247897382661923, 1.2722938635067185, 1.3203960685007963, 1.3690223971460211, 1.4181085678582204, 1.4675995755863493, 1.5174485540337082, 1.5676155237605893, 1.618066252352825, 1.6687712866664093, 1.7197051544565682], [0.3671818593760874, 0.4273099555248045, 0.48030520279500843, 0.5299157514716203, 0.5775026690234416, 0.6237348548957881, 0.6689956768144818, 0.7135268149734143, 0.7574910261402733, 0.8010039005559892, 0.8441531886668873, 0.8870159805810948, 0.9296819454693137, 0.9722847875737038, 1.0150257137074803, 1.0581560252499371, 1.1019108444708139, 1.14644132071577, 1.191795827318961, 1.237946348441428, 1.2848258998189876, 1.332355597388796, 1.3804587294116382, 1.4290662987793543, 1.47811829977904, 1.5275632294876373, 1.5773570274858204, 1.6274619429079007, 1.6778455089538895, 1.728479671620392], [0.3809039601281168, 0.441625357508203, 0.4950579772071751, 0.5450215293169391, 0.592905628384904, 0.6393943960608762, 0.684880579209215, 0.7296121558754571, 0.7737562922141495, 0.8174314047100478, 0.8607255084279288, 0.9037097200581007, 0.9464538310086525, 0.9890497672358259, 1.031640595238047, 1.074434085967645, 1.1176723674036657, 1.1615641532137952, 1.206230275727178, 1.251698172950591, 1.29793120131994, 1.3448621452326543, 1.392415753485645, 1.4405200107265037, 1.4891104117996206, 1.5381308129822053, 1.5875328849640258, 1.6372751250674609, 1.6873218283515787, 1.7376421613694726], [0.3946030664982533, 0.45590479317563604, 0.5097665401936732, 0.5600773534282057, 0.6082544853491102, 0.6549967906085674, 0.7007060955734936, 0.745636485770827, 0.7899594622938103, 0.8337965130063513, 0.8772372697350256, 0.9203513562803687, 0.9631993601813139, 1.005848396960646, 1.048395905470431, 1.0909948847219129, 1.1338567781986522, 1.177210856951684, 1.2212390801614885, 1.2660352462966569, 1.3116097425106372, 1.3579196955714976, 1.4048985022340292, 1.4524745125147644, 1.500580048592784, 1.5491546860828214, 1.5981457818370657, 1.6475078922219135, 1.6972018529768376, 1.7471938418714141], [0.40828060615170936, 0.4701503878417938, 0.5244333923506403, 0.5750859371449412, 0.6235520687583429, 0.6705449202924342, 0.7164751179100428, 0.7616026706748593, 0.8061033042038006, 0.8501016039229141, 0.8936893114443801, 0.9369365662653522, 0.9799005464704731, 1.0226356598795001, 1.065209851259014, 1.1077282616340083, 1.1503532866279824, 1.193296705301649, 1.2367719276696372, 1.2809346821731726, 1.325855062875235, 1.3715297446392003, 1.4179115172874004, 1.4649349903305449, 1.5125320291134454, 1.560638961860669, 1.609199081596174, 1.6581629285887822, 1.7074876913405472, 1.7571363493937837], [0.42193786930169463, 0.48436407116397207, 0.539060813233394, 0.5900497628044002, 0.6388009750574487, 0.6860414380854563, 0.7321903166172703, 0.7775133694144438, 0.8221904298869633, 0.8663491318422916, 0.9100834694056075, 0.953464912915314, 0.9965500285662346, 1.0393876629881849, 1.082029343600645, 1.124546555595792, 1.1670534218350024, 1.2097203227650335, 1.2527559220488846, 1.2963555518573164, 1.3406499110712036, 1.3856886169926448, 1.4314574627627306, 1.4779063849073477, 1.524971622838554, 1.5725884363463942, 1.6206968547560745, 1.6692435662677385, 1.7181820149790783, 1.7674717973091016], [0.4355760263958101, 0.4985476012761773, 0.5536508874824655, 0.6049711080649458, 0.6540035939206876, 0.7014887925579254, 0.7478541629656072, 0.793371050132741, 0.8382232815392785, 0.8825414677279481, 0.9264218658943423, 0.969937569789642, 1.0131457609624241, 1.056093429714234, 1.098824126067229, 1.141389125508906, 1.1838655307729158, 1.2263769184651785, 1.269099588087362, 1.3122357903105406, 1.3559613660249685, 1.400383670303986, 1.44553487000199, 1.4913922491547538, 1.5379040456757875, 1.5850084083926428, 1.6326438493050486, 1.6807538190672793, 1.7292881143413765, 1.778202835553653], [0.44919614288276927, 0.5127025850719835, 0.5682055271589075, 0.6198520686244691, 0.6691621303738328, 0.7168892491310782, 0.763468949034699, 0.8091780077153854, 0.8542041384718012, 0.8986808527034238, 0.9427066354885932, 0.9863562848893795, 1.0296880880487078, 1.0727489170280564, 1.1155790435923199, 1.15821909673683, 1.2007222597919356, 1.2431727700608055, 1.2857034088519084, 1.3284938929819807, 1.3717371165300702, 1.415588768367395, 1.4601348543567876, 1.50539285640928, 1.5513335187950765, 1.5979042724268655, 1.6450453520913038, 1.6926983680184804, 1.7408099359497209, 1.7893327051027785], [0.4627991917628816, 0.5268304956949847, 0.5827264908781407, 0.6346945777238131, 0.6842786240634932, 0.7322449084941749, 0.7790368052096889, 0.8249363799054301, 0.8701351293385163, 0.9147693914833401, 0.9589398302704358, 1.0027229492686447, 1.0461783133600595, 1.0893534243424712, 1.132287623749254, 1.1750166262684443, 1.2175800960070364, 1.2600349166640894, 1.3024734492527046, 1.34503680195391, 1.387905380150118, 1.4312602659076283, 1.4752369940199157, 1.5199023756357442, 1.5652616837202618, 1.6112807433707115, 1.6579068606986052, 1.705082456360253, 1.7527520798115654, 1.8008652774237743]]
zeros = [[0.05302558037376291, 0.10766084030247812, 0.16202764889502647, 0.21630297949929894, 0.27053710602765524, 0.3247493649037821, 0.3789486919989324, 0.433139760881445, 0.4873252445761774, 0.5415067790948483, 0.5956854204062967, 0.6498618798978972, 0.704036653996051, 0.7582100994844437, 0.8123824792847765, 0.866553991364633, 0.9207247875974955, 0.9748949864210336, 1.0290646815463889, 1.0832339480831203, 1.1374028469314015, 1.1915714279870633, 1.2457397325173694, 1.2999077949472895, 1.3540756442200172, 1.408243304845511, 1.4624107977173961, 1.5165781407557681, 1.570745349417675, 1.624912437105988], [0.0574091383391586, 0.11031590952427449, 0.16389552081180578, 0.21773510202147417, 0.27169551328186103, 0.32572077209022937, 0.37978456739568284, 0.43387303245343317, 0.4879782048622205, 0.5420952003605737, 0.596220862092304, 0.650353061473167, 0.7044903122361735, 0.7586315457508614, 0.8127759742825139, 0.8669230045707832, 0.9210721814098505, 0.9752231497624905, 1.0293756286836009, 1.08352939297612, 1.1376842600334274, 1.1918400802362599, 1.2459967298332315, 1.300154105587226, 1.3543121206973883, 1.4084707016559408, 1.4626297857992039, 1.5167893193804092, 1.57094925603914, 1.625109555575374], [0.06841073579953012, 0.11790099943696784, 0.16940009291883742, 0.2219962772068173, 0.27515516705443965, 0.3286269569352178, 0.3822875566352635, 0.43606994238108193, 0.4899351487852237, 0.5438591083411367, 0.5978262007759684, 0.6518258658145551, 0.7058507167120128, 0.7598954357768036, 0.81395609958288, 0.8680297513372389, 0.9221141211506383, 0.9762074379228901, 1.0303082997202186, 1.084415582509574, 1.1385283746514505, 1.192645929062775, 1.2467676277355453, 1.3008929550460335, 1.355021477417181, 1.4091528276389578, 1.4632866926489512, 1.5174228039146507, 1.5715609297938136, 1.62570086941435], [0.08247156971922939, 0.12938618237372707, 0.17824504002773836, 0.22898109555216878, 0.2808700018130878, 0.3334443407772762, 0.38644410835312804, 0.4397220086384576, 0.4931904158158231, 0.5467945343843704, 0.6004985436872038, 0.6542781173324181, 0.7081161886562366, 0.7620004461235597, 0.8159217927912684, 0.8698733654809147, 0.9238498910498554, 0.9778472524287642, 1.0318621890394748, 1.0858920855598921, 1.1399348201353248, 1.1939886534328916, 1.2480521462909266, 1.302124097733042, 1.3562034977125916, 1.4102894906645058, 1.4643813470902218, 1.518478441185582, 1.572580233065191, 1.626686254518961], [0.09744025929851802, 0.14340052866358116, 0.18993105061080637, 0.23850882659478947, 0.2887635325447038, 0.3401348363189222, 0.3922326559684076, 0.4448158077672976, 0.49773505039647997, 0.5508951848711029, 0.6042332873346686, 0.6577063405565672, 0.7112840361327591, 0.7649444471275252, 0.8186713388022155, 0.8724524446030076, 0.9262783292381673, 0.9801416195761865, 1.03403647233095, 1.0879581979850863, 1.1419029901321158, 1.1958677274029414, 1.2498498262966724, 1.3038471303141907, 1.3578578253802591, 1.411880374571035, 1.4659134672025211, 1.5199559787293828, 1.5740069388711364, 1.6280655060638545], [0.112486381331199, 0.15867155431177268, 0.2037623342173625, 0.25030259324061715, 0.29872305603931054, 0.3486455538568076, 0.399624451664856, 0.451333806537418, 0.5035574337359848, 0.5561528896880807, 0.6090244313170637, 0.6621059807525097, 0.7153507110131715, 0.7687246162472268, 0.822202452903388, 0.8757651119367899, 0.9293978744337814, 0.9830892261451708, 1.036830034251375, 1.0906129642950895, 1.144432059744404, 1.198282433839144, 1.252160040351474, 1.3060615007208447, 1.3599839720689277, 1.4139250452697678, 1.4678826653988555, 1.5218550690418238, 1.5758407344432548, 1.6298383415333402], [0.12739196654636548, 0.17435500508559035, 0.21894318892149162, 0.2639595975111968, 0.3105813993354397, 0.3589038190996785, 0.40858309437205803, 0.4592549517848324, 0.5106439610084527, 0.5625581522667766, 0.6148649874644674, 0.6674717008665381, 0.720312024086906, 0.7733375950408118, 0.8265123966138705, 0.8798091029028525, 0.9332066314437921, 0.9866884695846292, 1.0402415073296094, 1.0938552081477682, 1.1475210096835284, 1.201231883732758, 1.2549820084409804, 1.3087665208315664, 1.3625813276791698, 1.4164229593352347, 1.47028845557014, 1.5241752755578204, 1.5780812262613526, 1.6320044049838378], [0.14213401063875716, 0.190049356149086, 0.234773393953764, 0.2789586832827678, 0.3240872861429875, 0.3708032137424484, 0.4190605855166004, 0.46855418237245505, 0.5189794563236195, 0.5701006868057249, 0.6217474407971294, 0.6737977379842249, 0.7261634118461553, 0.7787796874196107, 0.8315981255377015, 0.8845818762889761, 0.9377024556243989, 0.9909375228735078, 1.0442693223181598, 1.0976835718477653, 1.1511686575650637, 1.204715041230001, 1.2583148181339328, 1.3119613829031374, 1.365649173834925, 1.4193734751305087, 1.4731302623351499, 1.5269160803910955, 1.5807279465680044, 1.6345632725579484], [0.15673213367741484, 0.20562006822246132, 0.2507943960149334, 0.29475145249236023, 0.33888774549362566, 0.3841795229415888, 0.4309869416661719, 0.4791993959979516, 0.5285467431288462, 0.5787697221458775, 0.6296641786876424, 0.681078287293665, 0.7329002427260486, 0.7850470946861076, 0.8374564677809696, 0.8900807497794218, 0.9428830564507779, 0.9958344143338256, 1.0489117702753958, 1.1020965651002226, 1.1553736965566903, 1.2087307545469257, 1.2621574495247776, 1.3156451798587143, 1.3691867005123513, 1.4227758665349903, 1.4764074324500047, 1.530076893877072, 1.5837803613978183, 1.6375144592780295], [0.1712090105131512, 0.2210429176302298, 0.26678311245929404, 0.31089288903469736, 0.35456203131466063, 0.3987878201180169, 0.44425159691947347, 0.49114373217290774, 0.5393243352285461, 0.5885536534287767, 0.6386077261588408, 0.689307851265063, 0.7405181389432854, 0.7921361784351519, 0.8440843298584747, 0.8963030590301578, 0.9487461202369529, 1.0013771227100368, 1.0541670767445843, 1.107092623356659, 1.1601347416636016, 1.2132777930011562, 1.2665088051078959, 1.3198169295797195, 1.373193025941784, 1.4266293393662386, 1.480119248449688, 1.533657065976661, 1.5872378801674623, 1.64085742715625], [0.18558334697099543, 0.23632679113764457, 0.2826603024197176, 0.32710792137149275, 0.37071057880094227, 0.4143078019416735, 0.4586816072492326, 0.5043113653876966, 0.5512806700213311, 0.5994382894808725, 0.6485704781314076, 0.6984814323233846, 0.7490132666630884, 0.8000437332604239, 0.8514789235135805, 0.9032463382127742, 0.9552894516485605, 1.0075636877661547, 1.0600334884499054, 1.112670176175441, 1.1654503840207076, 1.218354890465015, 1.2713677448278782, 1.3244756033951737, 1.3776672199355415, 1.430933050616239, 1.4842649446166425, 1.5376558996162366, 1.5910998668863559, 1.6445915946690097], [0.19986984543398115, 0.2514868208183108, 0.2984099754609048, 0.34326626183296727, 0.3870380840882659, 0.43039608117791545, 0.47403379609198404, 0.518578694398244, 0.5643632591332378, 0.6114025100870836, 0.659543369529148, 0.7085943361030427, 0.758382501437552, 0.8087672328223924, 0.8596379995894752, 0.9109085171325703, 0.9625111326298206, 1.014392336344345, 1.0665093728953325, 1.1188277260795951, 1.1713192536446446, 1.2239607956251655, 1.2767331266192483, 1.3296201590183094, 1.3826083308389137, 1.4356861306574624, 1.4888437253988414, 1.5420726660604662, 1.5953656530646758, 1.6487163476567783], [0.21408012747471158, 0.2665372735011808, 0.31403735202516364, 0.35932184677160406, 0.4033741270226148, 0.4467602180688064, 0.49002096344172036, 0.5337612430404778, 0.5784831030835891, 0.6244099479782005, 0.6715125798739913, 0.7196411622170635, 0.7686232906009226, 0.8183049771126465, 0.8685600604904888, 0.9192881241556878, 0.9704096950745803, 1.0218616229983002, 1.0735933310095618, 1.125563938378525, 1.1777400911603446, 1.23009432951389, 1.2826038528249095, 1.3352495782450058, 1.3880154163542893, 1.44088770861596, 1.4938547864326523, 1.5469066224404655, 1.600034552415512, 1.6532310517262392], [0.22822352762408288, 0.28149031086528753, 0.3295526831165529, 0.37526579652363695, 0.41964141144167666, 0.4632030943900126, 0.5063683738335687, 0.5496230969848249, 0.5935003098937532, 0.6383963197561063, 0.6844531217083019, 0.7316132864194841, 0.7797328878705956, 0.8286560025920752, 0.878244494616588, 0.9283844724856954, 0.9789842988040033, 1.029970582436875, 1.0812843222995137, 1.132877742226777, 1.1847118290126406, 1.2367544508377808, 1.2889789229577688, 1.3413629097992603, 1.3938875785486038, 1.4465369411577595, 1.4992973383649173, 1.5521570315905635, 1.6051058774722418, 1.658135066250656], [0.24230766606901794, 0.29635613307188513, 0.34496628582331623, 0.39110224627400775, 0.4358136430680817, 0.4796184515265025, 0.5228653576268905, 0.5659154993866123, 0.6092235338361354, 0.6532556351841163, 0.6983186346016776, 0.7444938906006104, 0.7917064184422128, 0.8398195029414177, 0.8886915242943252, 0.9381977849446799, 0.9882348969783948, 1.038718887118571, 1.0895818004277062, 1.1407684426321132, 1.192233682518514, 1.2439403296335276, 1.2958574933283644, 1.3479593177809703, 1.4002240034197448, 1.4526330449906344, 1.5051706337141957, 1.5578231843833779, 1.6105789582709662, 1.663427760087049], [0.2563388542288129, 0.3111433337424341, 0.3602873782318413, 0.4068387085985122, 0.45188670080365895, 0.4959604829607897, 0.5393814943890587, 0.582425514355524, 0.6254335535043775, 0.6688341606871996, 0.7130289765001084, 0.7582497729409033, 0.8045330056283775, 0.8517933352327219, 0.89990176926706, 0.9487291714005994, 0.9981623531382636, 1.0481069956781721, 1.0984858539307587, 1.1492358418451616, 1.2003052506518925, 1.2516514296563896, 1.3032389440709846, 1.3550381362268216, 1.4070240054625214, 1.4591753334437285, 1.511473997062645, 1.56390442479504, 1.6164531632816594, 1.6691085291555854], [0.2703223860649906, 0.3258592180426681, 0.3755239839783129, 0.42248290050114157, 0.467864328932611, 0.5122146579478928, 0.5558508819840924, 0.5990048402839252, 0.6419239255202814, 0.6849426397430846, 0.7284615714882243, 0.772820565149813, 0.818189206262092, 0.8645709863975012, 0.9118750894499345, 0.9599803021733302, 1.0087684413797977, 1.0581362634054645, 1.1079973405448782, 1.1582803659259258, 1.208926625390688, 1.259887599522933, 1.311122954001671, 1.3625989303089139, 1.4142870777273153, 1.4661632575476586, 1.5182068589257434, 1.5704001779799646, 1.6227279228106652, 1.6751768160579499], [0.28426275152342284, 0.3405100530537537, 0.39068307636997684, 0.4380418994073031, 0.48375237778749774, 0.5283797527553804, 0.5722467090643586, 0.6155684432847831, 0.6585352463179229, 0.7013872871724051, 0.744455651819139, 0.7881091904040054, 0.8326298319588319, 0.878136309267031, 0.9246081988783991, 0.971952509603131, 1.0200556047163498, 1.0688089597576318, 1.1181179918720978, 1.1679031871554242, 1.2180985062733658, 1.268649171745897, 1.319509583423221, 1.3706415646171428, 1.4220129478814338, 1.473596452087242, 1.5253687937026574, 1.577309981685926, 1.6294027551410473, 1.6816321319475045], [0.29816379596295034, 0.35510126090875993, 0.40577074556236725, 0.45352201644335977, 0.49955683647566596, 0.544459253535199, 0.5885614056548243, 0.6320759314286262, 0.6751656587666272, 0.7180040975032664, 0.7608339380319452, 0.8039808598194368, 0.8477785058025659, 0.892455807935068, 0.938090439566475, 0.9846449044092604, 1.0320262570138559, 1.0801280925832275, 1.1288504439281648, 1.178106323391495, 1.2278223125892445, 1.277937066987008, 1.3283993642556244, 1.3791662787075587, 1.430201640721425, 1.4814747871181178, 1.5329595621567165, 1.5846335213916356, 1.6364772957246432, 1.6884740809045162], [0.31202884132573316, 0.36963756860648145, 0.42079234310575836, 0.468928856058763, 0.5152832741298559, 0.5604579003863743, 0.6047954003267386, 0.6485119827898993, 0.6917606579706093, 0.73467768970622, 0.777433518337922, 0.8202767581292, 0.8635234594173152, 0.9074698181857087, 0.9522973490033833, 0.998050977285565, 1.044681291555428, 1.0920968661037715, 1.140198112540024, 1.188892678252304, 1.238100277783333, 1.2877528973613892, 1.3377933953615295, 1.3881737694838558, 1.4388535473661668, 1.48979842539075, 1.5409791588962083, 1.5923706695929705, 1.6439513297873352, 1.6957023871134405], [0.3258607796675965, 0.3841231255336655, 0.4357525998115927, 0.4842674082973177, 0.5309367381819143, 0.5763804835202234, 0.6209519164299439, 0.6648731001941482, 0.7082959114166948, 0.7513397606714731, 0.7941293620784216, 0.8368396743496629, 0.8797248725135823, 0.9230865554154205, 0.9671827067458888, 1.0121532623039264, 1.0580173434797007, 1.1047174973694465, 1.1521647734867886, 1.2002659563252094, 1.248935496599165, 1.2980990561696344, 1.3476934380382164, 1.3976652787195507, 1.4479695008463687, 1.498567885944632, 1.5494278652894589, 1.6005215296858692, 1.6518248287504143, 1.7033169251834233], [0.33966214635921843, 0.39856159674397923, 0.4506557210844976, 0.4995421352327843, 0.546521781348141, 0.5922314888741371, 0.63703489240775, 0.6811605860595553, 0.7247630179274056, 0.7679560168779361, 0.810840496381791, 0.8535389733876662, 0.896234666455255, 0.9391847213734926, 0.9826718143810313, 1.0269163257780487, 1.0720223760981362, 1.117989007268895, 1.1647536218970285, 1.2122303393628833, 1.2603318712665734, 1.308978769900466, 1.3581020013092915, 1.4076426816155028, 1.4575508567522208, 1.507784113705334, 1.5583063080991084, 1.6090864848672444, 1.6600979908935647, 1.711317753992176], [0.35343517809310554, 0.412956237813873, 0.46550546440334345, 0.5147570446287644, 0.562042516122013, 0.6080150284431043, 0.6530482503360591, 0.6973773819158529, 0.7411608196658254, 0.7845121534949302, 0.8275221488583415, 0.8702828297961274, 0.9129204767036966, 0.9556276789619204, 0.9986607732910802, 1.0422799117074963, 1.0866695948035818, 1.1319035768835655, 1.1779654870137382, 1.2247897382657538, 1.2722938635064929, 1.3203960685006908, 1.369022397145672, 1.4181085678578993, 1.4675995755860596, 1.5174485540334393, 1.5676155237603382, 1.6180662523526481, 1.6687712866662294, 1.7197051544564603], [0.36718185936876874, 0.42730995552263945, 0.48030520278784394, 0.5299157514643009, 0.5775026690156936, 0.623734854891891, 0.6689956768082941, 0.7135268149695712, 0.7574910261396898, 0.8010039005541509, 0.8441531886640962, 0.8870159805807907, 0.9296819454686077, 0.9722847875717988, 1.0150257137056131, 1.058156025249762, 1.1019108444701011, 1.1464413207148638, 1.1917958273184972, 1.2379463484408364, 1.2848258998187934, 1.332355597388356, 1.3804587294113242, 1.4290662987790594, 1.4781182997787745, 1.527563229487434, 1.5773570274857922, 1.6274619429078476, 1.6778455089536608, 1.7284796716202027], [0.380903960125785, 0.4416253574963793, 0.49505797720617384, 0.545021529307, 0.5929056283765574, 0.6393943960531718, 0.684880579202118, 0.7296121558719894, 0.7737562922092197, 0.8174314047048697, 0.8607255084231129, 0.9037097200546526, 0.9464538310065533, 0.989049767233548, 1.0316405952355068, 1.0744340859669836, 1.117672367402148, 1.1615641532131458, 1.206230275726407, 1.251698172950156, 1.2979312013195443, 1.3448621452323999, 1.39241575348532, 1.4405200107264877, 1.489110411799292, 1.538130812982019, 1.587532884963912, 1.6372751250673485, 1.687321828351443, 1.7376421613693518], [0.3946030664944521, 0.4559047931671091, 0.5097665401823587, 0.5600773534189374, 0.6082544853400039, 0.6549967906031109, 0.7007060955709087, 0.7456364857640391, 0.789959462287673, 0.8337965130007573, 0.8772372697310418, 0.9203513562762087, 0.9631993601776074, 1.0058483969572924, 1.048395905469413, 1.0909948847209165, 1.1338567781973126, 1.1772108569509463, 1.2212390801611472, 1.2660352462959894, 1.3116097425100914, 1.3579196955710042, 1.4048985022335407, 1.452474512514345, 1.5005800485927159, 1.5491546860825391, 1.598145781836863, 1.6475078922218092, 1.697201852976712, 1.7471938418712902], [0.4082806061370169, 0.4701503878293404, 0.5244333923396444, 0.5750859371425149, 0.6235520687559654, 0.6705449202899638, 0.7164751179068549, 0.7616026706686864, 0.8061033041984039, 0.8501016039173565, 0.8936893114398526, 0.9369365662606642, 0.9799005464662317, 1.022635659876101, 1.0652098512573505, 1.107728261631866, 1.1503532866261303, 1.1932967053001597, 1.2367719276692641, 1.2809346821722414, 1.325855062875037, 1.371529744638689, 1.4179115172868484, 1.4649349903305269, 1.5125320291134, 1.5606389618606171, 1.6091990815960788, 1.6581629285887054, 1.7074876913403225, 1.7571363493935623], [0.421937869294669, 0.48436407116047514, 0.5390608132263066, 0.5900497627968386, 0.63880097505653, 0.6860414380769997, 0.7321903166104856, 0.7775133694077616, 0.8221904298805269, 0.8663491318395411, 0.9100834694001417, 0.9534649129137267, 0.9965500285657335, 1.0393876629845709, 1.082029343597425, 1.1245465555927474, 1.1670534218324788, 1.209720322763197, 1.2527559220484084, 1.2963555518560443, 1.3406499110708798, 1.3856886169919012, 1.4314574627620993, 1.47790638490685, 1.5249716228381347, 1.5725884363460132, 1.6206968547559042, 1.6692435662674372, 1.7181820149790616, 1.767471797308936], [0.43557602639421394, 0.4985476012633913, 0.5536508874777544, 0.6049711080608624, 0.654003593911428, 0.7014887925521709, 0.7478541629612375, 0.7933710501313415, 0.8382232815338192, 0.8825414677217777, 0.9264218658917028, 0.9699375697845453, 1.0131457609589805, 1.0560934297100393, 1.0988241260655423, 1.1413891255074278, 1.1838655307700332, 1.2263769184644555, 1.2690995880854725, 1.3122357903095232, 1.3559613660238605, 1.400383670303163, 1.4455348700016557, 1.491392249154304, 1.5379040456757003, 1.5850084083922293, 1.6326438493048576, 1.6807538190670834, 1.7292881143412557, 1.778202835553504], [0.44919614287513066, 0.5127025850589377, 0.5682055271470254, 0.6198520686217264, 0.6691621303694132, 0.7168892491244149, 0.7634689490331181, 0.8091780077151554, 0.8542041384684989, 0.898680852700221, 0.9427066354831555, 0.9863562848849856, 1.0296880880471224, 1.0727489170266877, 1.1155790435916426, 1.1582190967355628, 1.200722259789491, 1.2431727700588815, 1.2857034088497306, 1.3284938929812884, 1.3717371165294838, 1.4155887683666046, 1.460134854356354, 1.5053928564089374, 1.5513335187945074, 1.5979042724264862, 1.6450453520909294, 1.692698368018146, 1.7408099359496436, 1.789332705102538], [0.46279919175317596, 0.5268304956814541, 0.5827264908661104, 0.6346945777132167, 0.6842786240542162, 0.7322449084911644, 0.7790368052044662, 0.8249363798981983, 0.8701351293353243, 0.9147693914771019, 0.958939830267028, 1.0027229492675704, 1.046178313355645, 1.0893534243378793, 1.1322876237452366, 1.1750166262648287, 1.2175800960058427, 1.2600349166631726, 1.3024734492501804, 1.3450368019525742, 1.3879053801485304, 1.4312602659065952, 1.4752369940198908, 1.5199023756349779, 1.5652616837196929, 1.611280743370298, 1.6579068606983858, 1.7050824563598612, 1.7527520798114635, 1.8008652774235256]]

number_zeros = 30

def Beta(m, n):
    return zeros[m][n - 1]

def R(r, m, n):
    return sp.yv(m,a*Beta(m, n))*sp.jv(m,r*Beta(m, n)) - sp.jv(m,a*Beta(m, n))*sp.yv(m,r*Beta(m, n))

def N(m, n):
    return (2/((((np.pi)**2))*((Beta(m, n))**2)))*((((sp.jv(m, a*Beta(m, n)))**2)/(((sp.jv(m, b*Beta(m, n)))**2))) - 1)

def R_total(r, r_prime, m, n):
    return (R(r, m, n)*R(r_prime, m, n))/N(m, n)

def Delta_function(r, r_prime, m):
    count = 1
    nascent_delta = 0
    while count <= number_zeros:
        nascent_delta += R_total(r, r_prime, m, count)
        count += 1
    return nascent_delta

r = a

r_list_0 = []
delta_list_0 = []

r_list_5 = []
delta_list_5 = []

r_list_10 = []
delta_list_10 = []

r_list_15 = []
delta_list_15 = []

r_list_16 = []
delta_list_16 = []

r_list_20 = []
delta_list_20 = []

r_list_23 = []
delta_list_23 = []

r_list_24 = []
delta_list_24 = []

r_prime = 30

count = 1

while r <= b:
    r_list_0.append(r)
    r_list_5.append(r)
    r_list_10.append(r)
    r_list_15.append(r)
    r_list_16.append(r)
    r_list_20.append(r)
    r_list_23.append(r)
    r_list_24.append(r)
    
    delta_0 = Delta_function(r, r_prime, 0)
    delta_list_0.append(delta_0)

    delta_5 = Delta_function(r, r_prime, 5)
    delta_list_5.append(delta_5)

    delta_10 = Delta_function(r, r_prime, 10)
    delta_list_10.append(delta_10)

    delta_15 = Delta_function(r, r_prime, 15)
    delta_list_15.append(delta_15)

    delta_16 = Delta_function(r, r_prime, 16)
    delta_list_16.append(delta_16)

    delta_20 = Delta_function(r, r_prime, 20)
    delta_list_20.append(delta_20)

    delta_23 = Delta_function(r, r_prime, 23)
    delta_list_23.append(delta_23)

    delta_24 = Delta_function(r, r_prime, 24)
    delta_list_24.append(delta_24)

    r += 1
    print(count)
    count += 1

plt.title("Delta Functions at Ten Digits of Precision")
plt.xlabel("r (cm)")
plt.ylabel("Delta Function")

plt.plot(r_list_0, delta_list_0, label = "0")
plt.plot(r_list_5, delta_list_5, label = "5")
plt.plot(r_list_10, delta_list_10, label = "10")
plt.plot(r_list_15, delta_list_15, label = "15")
plt.plot(r_list_16, delta_list_16, label = "16")
plt.plot(r_list_20, delta_list_20, label = "20")
plt.plot(r_list_23, delta_list_23, label = "23")
#plt.plot(r_list_24, delta_list_24, label = "25")

plt.legend()
plt.show()
