
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;

//To draw the scaffold map with IR args[2]
public class getScafA {
	public static void main(String[] args) throws IOException, IOException {
		
		//args[0] is the path for $iA.pos files
		//args[1] is the blocks.txt file with its path
		//args[2] is the representative letter for IR
		//args[3] is the total number of genomes
		//args[4] is the length of E
		//args[5] is the threshold for short E
		
		
		String opName=args[0]+"/scafMapWith"+args[2]+".txt";
		PrintWriter writer=new PrintWriter(opName,"UTF-8"); //like the mgr_macro.txt file
		PrintWriter writer1=new PrintWriter(args[0]+"/blks_"+args[2]+".txt","UTF-8"); // like the blocks.txt file
	    writer1.println("#scaffold_rep start length strand startPos_on_A");	
		ArrayList<scaffold> scafs;
		for(int i=0; i<Integer.parseInt(args[3]);i++){
			writer.println(">Strain "+(i+1));
			writer1.println(">Strain "+(i+1));
			scafs=new ArrayList<scaffold>(); //To store all the scaffolds for strain i
			
			File posf = new File(args[0]+"/"+(i+1)+args[2]+".pos");
			FileReader fr = new FileReader(posf);  
			BufferedReader br = new BufferedReader(fr);
			String line;
			String[] tokens;
			int s;
			int l;
			char strand;
			int refs;
			while(true){
				line=br.readLine();
				if(line==null) break;
				if(line.contains("No")) continue;
				tokens=line.split(" ");
				if(Integer.parseInt(tokens[0])<=Integer.parseInt(tokens[1])){
					s=Integer.parseInt(tokens[0]);
					l=Integer.parseInt(tokens[1])-Integer.parseInt(tokens[0])+1;
					strand='+';
				}
				else{
					s=Integer.parseInt(tokens[1]);
					l=Integer.parseInt(tokens[0])-Integer.parseInt(tokens[1])+1;
					strand='-';
				}
				refs=Integer.parseInt(tokens[2]);
				if(l<Double.parseDouble(args[5])*Integer.parseInt(args[4]))
					scafs.add(new scaffold(args[2]+"s",s,l,strand,refs));
				else
					scafs.add(new scaffold(args[2],s,l,strand,refs));
			} //Finish adding the IR's position in strain i+1
			
			File blocksf = new File(args[1]);
			FileReader fr1 = new FileReader(blocksf);  
			BufferedReader br1 = new BufferedReader(fr1);
	
			while(true){
				line=br1.readLine();
				if(line==null) break;
				if(line.contains("#")) continue;
				tokens=line.split(" ");
				s=Integer.parseInt(tokens[4*i+2]);
				l=Integer.parseInt(tokens[4*i+3]);
				if(tokens[4*i+4].equals("1"))
					strand='+';
				else
					strand='-';
				scafs.add(new scaffold(tokens[0],s,l,strand,0));
			}
			
			scafs.sort(new CustomComparator());
			
			for(int k=0;k<scafs.size();k++){
				writer.print(scafs.get(k).strand+scafs.get(k).rep+" ");
                if(k==scafs.size()-1) writer.println("");
				writer1.println(scafs.get(k).rep+" "+scafs.get(k).start+" "+scafs.get(k).length+" "+scafs.get(k).strand+" "+scafs.get(k).refs);
			}
			fr.close();
			br.close();
			fr1.close();
			br1.close();
		}
		writer.close();
		writer1.close();
	}

}

//both core sequence and IR sequence are scaffold 
class scaffold{ 
	String rep; // for core sequence, rep is the id in mgr_macro.txt; For IR, rep is A/B/C/D...
	int start;
	int length;
	char strand;
	int refs; // for IR refs is the start pos in IR itself. for core seq, refs=0;
	//If two IR's refs and length are almost the same, they can be regarded as the same IR
	public scaffold(String rep,int s,int l,char strand,int refs){
		this.rep=rep;
		start=s;
		length=l;
		this.strand=strand;
		this.refs=refs;
	}
}

class CustomComparator implements Comparator<scaffold> {
    public int compare(scaffold o1, scaffold o2) {
        return o1.start - o2.start;
    }
}
