


public class getPrevBlk {
	
	//args[0] is the input blk which want to get the previous block
	//args[1-n] is the genome permutation
	public static void main(String[] args) {
		String target=args[0];
		int n=args.length;
		
		if(target.equals(args[1]))
			System.out.println("0"); // 0 means no previous block;
		else{
			for(int i=2;i<=n-2;i++){
				if(args[i].equals(target)){
					System.out.println(args[i-1]);
					break;
				}
			}
		}
	}

}
