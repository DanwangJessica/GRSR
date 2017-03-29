import java.util.Arrays;


public class getNbD {

		//args[0] is the start block of the rearrangement region on s
		//args[1] is the end block of the rearrangement region on s
		//args[2 to n-1] is the destination permutation d ended with $
		public static void main(String[] args) {
			String targetS=args[0];
			String targetE=args[1];
			int n=args.length;
			
			String[] pd=new String[n-3]; //pd is the permutation of d
			for(int i=0;i<pd.length;i++){
				pd[i]=args[i+2];
			}
			int dL=n-3; //size of permutation d
			
			int idS=Arrays.asList(pd).indexOf(targetS); //index for start block on d
			int idE=Arrays.asList(pd).indexOf(targetE);
			
			if(idS>=0 && idE>=0){ // targetS~ targetE on d 
				if(idS!=0 && idE!=dL-1){
					System.out.println(pd[idS-1]+" "+pd[idS]+" "+pd[idE]+" "+pd[idE+1]);
				}
				else{
					System.out.println("0 0 0 0");
				}
			}
			
			else if(idS<0 && idE<0){ //check -targetE ~ -targetS on d
				String revS=Integer.parseInt(targetS)*-1+"";
				String revE=Integer.parseInt(targetE)*-1+"";
				idS=Arrays.asList(pd).indexOf(revS); //index for start block on d
				idE=Arrays.asList(pd).indexOf(revE);
				if(idS>=0 && idE>=0){ 
					if(idE!=0 && idS!=dL-1){
						System.out.println(pd[idE-1]+" "+pd[idE]+" "+pd[idS]+" "+pd[idS+1]);
					}
					else{
						System.out.println("0 0 0 0");
					}
				}
				
			}
			
			else
				System.out.println("0 0 0 0");
			
		}

}
