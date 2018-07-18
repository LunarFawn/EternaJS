import {AppMode} from "../../flashbang/core/AppMode";
import {MissionClearedPanel} from "../mode/PoseEdit/MissionClearedPanel";
import {Background} from "../vfx/Background";

export class PoseTestMode extends AppMode {
    protected setup(): void {
        super.setup();

        this.addObject(new Background(), this.modeSprite);

        let panel = new MissionClearedPanel(false, PoseTestMode.MORE_INFO);
        this.addObject(panel, this.modeSprite);

        panel.createRankScroll(JSON.parse(PoseTestMode.SUBMISSION_RSP));
    }

    private static readonly SUBMISSION_RSP: string = '{"points-new":100,"pointsrank-before":{"richer":[{"name":"mayaa25","rank":108102,"points":"25"},{"name":"StumpyTheStump","rank":108102,"points":"25"},{"name":"adammay451","rank":108102,"points":"25"},{"name":"AsianSpartan","rank":108102,"points":"25"},{"name":"shub","rank":108102,"points":"25"},{"name":"IIDatAsheII","rank":108102,"points":"25"},{"name":"I\'maperson","rank":108102,"points":"25"},{"name":"wolfheart11","rank":108102,"points":"25"},{"name":"RedRavin","rank":108102,"points":"25"},{"name":"cambear","rank":108102,"points":"25"},{"name":"InderSharma4","rank":108102,"points":"25"},{"name":"RekerNAtaor5000","rank":108102,"points":"25"},{"name":"smirick","rank":108102,"points":"25"},{"name":"adamhanul","rank":108102,"points":"25"},{"name":"BlackElectrum","rank":108102,"points":"25"},{"name":"ScubaRW4","rank":108102,"points":"25"},{"name":"lanwill","rank":108102,"points":"25"},{"name":"Norman_RKN","rank":108102,"points":"25"},{"name":"jai289","rank":108102,"points":"25"},{"name":"LethalPianist","rank":108102,"points":"25"},{"name":"platoisawesome","rank":108102,"points":"25"},{"name":"bluewolf","rank":108102,"points":"25"},{"name":"Klykows","rank":108102,"points":"25"},{"name":"awoolf19","rank":108102,"points":"25"},{"name":"jmaxime","rank":108102,"points":"25"},{"name":"justin.parker@portageps.org","rank":108102,"points":"25"},{"name":"Christiaan","rank":108102,"points":"25"},{"name":"EbbeS","rank":108102,"points":"25"},{"name":"ninjakey","rank":108102,"points":"25"},{"name":"alibek.karimov","rank":108102,"points":"25"},{"name":"crystalrose_cx_","rank":108102,"points":"25"},{"name":"rebecca.ha","rank":108102,"points":"25"},{"name":"atnguyen02","rank":108102,"points":"25"},{"name":"krizzy","rank":108102,"points":"25"},{"name":"badriakafala","rank":108102,"points":"25"},{"name":"Sagemyster","rank":108102,"points":"25"},{"name":"Smeagol","rank":108102,"points":"25"},{"name":"Optkenate","rank":108102,"points":"25"},{"name":"adas16","rank":108102,"points":"25"},{"name":"mrfish","rank":108102,"points":"25"},{"name":"AustinMullins115","rank":108102,"points":"25"},{"name":"Kiongbrent","rank":108102,"points":"25"},{"name":"claybaby97","rank":108102,"points":"25"},{"name":"biologyclass2","rank":108102,"points":"25"},{"name":"EderPort","rank":108102,"points":"25"},{"name":"Dan Schuler","rank":108102,"points":"25"},{"name":"Rhubabr88","rank":108102,"points":"25"},{"name":"Cannon","rank":108102,"points":"25"},{"name":"nikkif","rank":108102,"points":"25"},{"name":"Kasey Anderson","rank":108102,"points":"25"}],"poorer":[{"name":"MildaZ","rank":108765,"points":null},{"name":"Mal.com","rank":108765,"points":null},{"name":"Thotknectar","rank":108765,"points":null},{"name":"max_dor","rank":108765,"points":null},{"name":"tho","rank":108765,"points":null},{"name":"silvaoliver","rank":108765,"points":null},{"name":"minijam3","rank":108765,"points":null},{"name":"megkelly","rank":108765,"points":null},{"name":"tybiousd","rank":108765,"points":null},{"name":"elaha.NAZ","rank":108765,"points":null},{"name":"Bernstein ","rank":108765,"points":null},{"name":"MurilloJose","rank":108765,"points":null},{"name":"apache","rank":108765,"points":null},{"name":"trevor.kaye","rank":108765,"points":null},{"name":"jewelryvang","rank":108765,"points":null},{"name":"loic621","rank":108765,"points":null},{"name":"Juuliana22","rank":108765,"points":null},{"name":"kaira","rank":108765,"points":null},{"name":"JoonBug","rank":108765,"points":null},{"name":"loneamigo","rank":108765,"points":null},{"name":"Jespez8","rank":108765,"points":null},{"name":"joaquinv","rank":108765,"points":null},{"name":"yayasoba","rank":108765,"points":null},{"name":"parker3286","rank":108765,"points":null},{"name":"allykinaman","rank":108765,"points":null},{"name":"crkdmn","rank":108765,"points":null},{"name":"bpethel","rank":108765,"points":null},{"name":"mtruo008","rank":108765,"points":null},{"name":"jwaldo","rank":108765,"points":null},{"name":"Leslye Esteves ","rank":108765,"points":null},{"name":"BrenBren","rank":108765,"points":null},{"name":"erinhaus","rank":108765,"points":null},{"name":"acep95","rank":108765,"points":null},{"name":"jorjohnson","rank":108765,"points":null},{"name":"shivi","rank":108765,"points":null},{"name":"cmkreiner","rank":108765,"points":null},{"name":"redesmond","rank":108765,"points":null},{"name":"farah_j","rank":108765,"points":null},{"name":"clsmith293","rank":108765,"points":null},{"name":"lizzy47","rank":108765,"points":null},{"name":"karthik.sivadas","rank":108765,"points":null},{"name":"quinn.botkin@puukumu.org","rank":108765,"points":null},{"name":"jeloriaga","rank":108765,"points":null},{"name":"Great_Guy96","rank":108765,"points":null},{"name":"imagahub","rank":108765,"points":null},{"name":"jul059","rank":108765,"points":null},{"name":"acio","rank":108765,"points":null},{"name":"jlamb23","rank":108765,"points":null},{"name":"wassup b#########","rank":108765,"points":null},{"name":"smiller23","rank":108765,"points":null}],"rank":108765,"points":0},"pointsrank-after":{"richer":[{"name":"champagne","rank":106785,"points":"110"},{"name":"doratheextruder","rank":106785,"points":"110"},{"name":"goda","rank":106785,"points":"110"},{"name":"silvano","rank":106785,"points":"110"},{"name":"HoxC6","rank":106785,"points":"110"},{"name":"snowflake85","rank":106785,"points":"110"},{"name":"noush\'","rank":106785,"points":"110"},{"name":"Houyo","rank":106785,"points":"110"},{"name":"oligogene","rank":106785,"points":"110"},{"name":"stevewelch","rank":106785,"points":"110"},{"name":"settler1000","rank":106785,"points":"110"},{"name":"cooleel","rank":106785,"points":"110"},{"name":"fa12","rank":106785,"points":"110"},{"name":"juanjomeister","rank":106785,"points":"110"},{"name":"lys350","rank":106785,"points":"110"},{"name":"blackbox","rank":106785,"points":"110"},{"name":"kal","rank":106785,"points":"110"},{"name":"synthesis laboratory","rank":106785,"points":"110"},{"name":"egzekutor12","rank":106785,"points":"110"},{"name":"18jarnold","rank":106785,"points":"110"},{"name":"hailee","rank":106785,"points":"110"},{"name":"AdymRNA","rank":106785,"points":"110"},{"name":"zackyshadow","rank":106785,"points":"110"},{"name":"stromd","rank":106785,"points":"110"},{"name":"iapurgill","rank":106785,"points":"110"},{"name":"laisar","rank":106785,"points":"110"},{"name":"raptor2910","rank":106785,"points":"110"},{"name":"jeamsbaugh","rank":106785,"points":"110"},{"name":"Courtney Lentz","rank":106785,"points":"110"},{"name":"rayrodriguez","rank":106785,"points":"110"},{"name":"Buckiguy","rank":106785,"points":"110"},{"name":"bcoleman15","rank":106785,"points":"110"},{"name":"magdi","rank":106785,"points":"110"},{"name":"keith.lundquist","rank":106785,"points":"110"},{"name":"username123","rank":106785,"points":"110"},{"name":"eliejdl","rank":106785,"points":"110"},{"name":"herbert","rank":106785,"points":"110"},{"name":"MikeHawk","rank":106785,"points":"110"},{"name":"feldyz","rank":106785,"points":"110"},{"name":"szeki2","rank":106785,"points":"110"},{"name":"porcbolsevic","rank":106785,"points":"110"},{"name":"andoresu","rank":106785,"points":"110"},{"name":"vero","rank":106785,"points":"110"},{"name":"ginger bret man","rank":106785,"points":"110"},{"name":"Madison.lee","rank":106785,"points":"110"},{"name":"cmtaylor01","rank":106961,"points":"108"},{"name":"samanthakqueved","rank":106962,"points":"105"},{"name":"jbenda","rank":106962,"points":"105"},{"name":"Flangan123","rank":106962,"points":"105"},{"name":"d_officeboy","rank":106962,"points":"105"}],"poorer":[{"name":"peters4oz_eterna","rank":106966,"points":"100"},{"name":"darrans","rank":106966,"points":"100"},{"name":"WillHines222","rank":106966,"points":"100"},{"name":"benS","rank":106966,"points":"100"},{"name":"Keaton Norton","rank":106966,"points":"100"},{"name":"Shen1990","rank":106966,"points":"100"},{"name":"z0mbiecraze","rank":106966,"points":"100"},{"name":"trinsa","rank":106966,"points":"100"},{"name":"fjhguifxh","rank":106966,"points":"100"},{"name":"bandlife","rank":106966,"points":"100"},{"name":"Marquis ","rank":106966,"points":"100"},{"name":"lt4191998","rank":106966,"points":"100"},{"name":"spagettigaling","rank":106966,"points":"100"},{"name":"dagr9029","rank":106966,"points":"100"},{"name":"Ctrudell5","rank":106966,"points":"100"},{"name":"mike0308","rank":106966,"points":"100"},{"name":"rickynaps","rank":106966,"points":"100"},{"name":"timesen","rank":106966,"points":"100"},{"name":"cyclopian","rank":106966,"points":"100"},{"name":"Aurora411","rank":106966,"points":"100"},{"name":"okc23","rank":106966,"points":"100"},{"name":"tsmith","rank":106966,"points":"100"},{"name":"foxtrot","rank":106966,"points":"100"},{"name":"dmainbeast040","rank":106966,"points":"100"},{"name":"DrGary","rank":106966,"points":"100"},{"name":"nicolebelli","rank":106966,"points":"100"},{"name":"isaiah3.0","rank":106966,"points":"100"},{"name":"pikachugo","rank":106966,"points":"100"},{"name":"Pianous","rank":106966,"points":"100"},{"name":"boomod","rank":106966,"points":"100"},{"name":"heilei","rank":106966,"points":"100"},{"name":"mjas21","rank":106966,"points":"100"},{"name":"kingofawe","rank":106966,"points":"100"},{"name":"quentin78700","rank":106966,"points":"100"},{"name":"Gordon.matthew.p","rank":106966,"points":"100"},{"name":"fcort48","rank":106966,"points":"100"},{"name":"thesortinghat13","rank":106966,"points":"100"},{"name":"WillebrandtC","rank":106966,"points":"100"},{"name":"jacodo","rank":106966,"points":"100"},{"name":"Turtle60","rank":106966,"points":"100"},{"name":"Ehawk726","rank":106966,"points":"100"},{"name":"kwontino","rank":106966,"points":"100"},{"name":"satpathya","rank":106966,"points":"100"},{"name":"jschuldinger","rank":106966,"points":"100"},{"name":"Lior Pyt","rank":106966,"points":"100"},{"name":"JillB","rank":106966,"points":"100"},{"name":"ctagoe","rank":106966,"points":"100"},{"name":"Deeps","rank":106966,"points":"100"},{"name":"abemonkeyman","rank":106966,"points":"100"},{"name":"SteamBlast","rank":106966,"points":"100"}],"rank":106966,"points":100},"next-puzzle":null,"success":true}';

    private static readonly MORE_INFO: string = "Long stacks and short stacks behave differently. Long stacks are generally more stable. Short stacks usually require more <font color='#B93F3C'>G</font>-<font color= '#4FB748'>C</font> pairs.<br><br><img width='274' hspace='55' height='174' src='/puzzle-progression/images/short-stack-long-stack.png'><br><br><br><br><br><br><br><br><br><br><br><br><br><p><font size='-3'><u><a target='_blank' href='http://www.eternagame.org/web/puzzle/6502997/'>Feedback</a></u></font><p>";
}

interface PoseDesc {
    seq: number[];
    barcodes?: number[];
    oligos?: number[];
    oligo?: number[];
    pairs: number[];
    structConstraints?: boolean[];
    puzlocks: boolean[];
    shiftLimit?: number;
    restrictedSequence?: number[];
}
