import Folder, {MultiFoldResult} from '../Folder';
import NuPACK from '../NuPACK';
//import Vienna from '../Vienna';
//import Vienna2 from '../Vienna2';
//import LinearFoldV from '../LinearFoldV';
import './jest-matcher-deep-close-to';
import SecStruct from 'eterna/rnatypes/SecStruct';
import Sequence from 'eterna/rnatypes/Sequence';
import {Oligo} from 'eterna/rnatypes/Oligo';

function CreateFolder(type: any): Promise<Folder | null> {
    return type.create();
}




test('NuPACK:multifold', () => {        
    return expect(CreateFolder(NuPACK)
      .then((folder) => {
       if (folder === null) return;
       

    //   public multifold(
    //    _seq: Sequence, _secstruct: SecStruct | null, _oligos: Oligo[],
    //    _desiredPairs: string | null = null, _temp: number = 37
    ////): MultiFoldResult | undefined {
    //    return undefined;
   //}


         //let sequence = Sequence.fromSequenceString("cgcaaccuac");
         let sequence = Sequence.fromSequenceString("UCGGAACUUAGCUUAGAUGGUUGCGUUGAAUUCGAGAUCUACAUGGUAGUUCGCUAUCAUGUAGAUUUCGGGUUCCAUCUGCAGU");
         //2,3,2,1,1,2,2,4,1,2
         //should represent CGCAACCAUC
         let oligos: Oligo[] = new Array(1);
         //oligos[0].sequence=Sequence.fromSequenceString("UCGGAACUUAGCUUAGAUGGUUGCGUUGAAUUCGAGAUCUACAUGGUAGUUCGCUAUCAUGUAGAUUUCGGGUUCCAUCUGCAGU");
         oligos[0] = {sequence: [2,3,2,1,1,2,2,4,1,2], malus:0};
         //oligos[0] = {sequence: [4, 2, 3, 3, 1, 1, 2, 4, 4, 1, 3, 2, 4, 4, 1, 3, 1, 4, 3, 3, 4, 4, 3, 2, 3, 4, 4, 3, 1, 1, 4, 4, 2, 3, 1,
         //                    3, 1, 4, 2, 4, 1, 2, 1, 4, 3, 3, 4, 1, 3, 4, 4, 2, 3, 2, 4, 1, 4, 2, 1, 4, 3, 4, 1, 3, 1, 4, 4, 4, 2, 3, 3, 3
         //                   , 4, 4, 2, 2, 1, 4, 2, 4, 3, 2, 1, 3, 4], malus: 0};
         //UCGGAACUUAGCUUAGAUGGUUGCGUUGAAUUCGAGAUCUACAUGGUAGUUCGCUAUCAUGUAGAUUUCGGGUUCCAUCUGCAGU
         //let temperature: number = 37;
           //let isPsuedoknot:boolean = false;
        let desiredPairs = null;
        let secstruct = null;
         
        
          const MultiFoldResultObject: 
          MultiFoldResult | undefined = folder.multifold (
            sequence, secstruct, oligos ,desiredPairs, 37  
           );	   
         
        
        
        let pairs: SecStruct = (MultiFoldResultObject as MultiFoldResult).pairs;
        
        console.log(pairs);
        
        expect(pairs).toBeDefined();
     

   }))
   .resolves.toBeUndefined(); // (we're returning a promise)

});