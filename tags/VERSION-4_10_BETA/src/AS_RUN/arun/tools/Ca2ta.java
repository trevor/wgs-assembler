import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.*;
import java.lang.reflect.Field;

class Compression {

	public static final byte[] compress(String str) throws IOException {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		ZipOutputStream zout = new ZipOutputStream(out);
		zout.putNextEntry(new ZipEntry("0"));
		zout.write(str.getBytes());
		zout.closeEntry();
		byte[] compressed = out.toByteArray();
		zout.close();
		return compressed;
	}

	public static final String decompress(byte[] compressed) throws IOException {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		ByteArrayInputStream in = new ByteArrayInputStream(compressed);
		ZipInputStream zin = new ZipInputStream(in);
		ZipEntry entry = zin.getNextEntry();
		byte[] buffer = new byte[1024];
		int offset = -1;
		while((offset = zin.read(buffer)) != -1) {
			out.write(buffer, 0, offset);
		}
		String decompressed = out.toString();
		out.close();
		zin.close();
		return decompressed;
	}

}

class ToString {
	/**
	* Intended only for debugging.
	*
	* <P>Here, a generic implementation uses reflection to print
	* names and values of all fields <em>declared in this class</em>. Note that
	* superclass fields are left out of this implementation.
	*
	* <p>The format of the presentation could be standardized by using
	* a MessageFormat object with a standard pattern.
	*/
	public String toString() {
	  StringBuilder result = new StringBuilder();

	  result.append( this.getClass().getName() );
	  result.append( "$ Object {" );
	  result.append("\n");

	  //determine fields declared in this class only (no fields of superclass)
	  Field[] fields = this.getClass().getDeclaredFields();

	  //print field names paired with their values
	  for ( Field field : fields  ) {
	    result.append("  ");
	    try {
              result.append( field.getName() );
              result.append(": ");
              //requires access to private field:
              result.append( field.get(this) );
	    }
	    catch ( IllegalAccessException ex ) {
              System.out.println(ex);
	    }
	    result.append("\n");
	  }
	  result.append("}");

	  return result.toString();
	}


}

class Contig {
	String acc;
	int len;
	int npc;
	byte [] cns;
	ArrayList<MPS> mpsList;
	
	public String toString() {
		StringBuilder result = new StringBuilder();
	     	try {
			String cnsDecompressed = Compression.decompress(cns);
			System.out.println("Compressed size: " + cns.length + ", Decompressed size: " + cnsDecompressed.getBytes().length);
			result.append("##" + acc + " " + npc + " " + cnsDecompressed.length() + " bases, 00000000 checksum." + "\n");
			result.append(Ca2ta.formatSeq(cnsDecompressed));
			Comparator<MPS> comp = new Comparator<MPS>() { public int compare(MPS o1, MPS o2) { return o1.mid.compareTo(o2.mid);}};
			Collections.sort(mpsList,comp);
			for ( MPS mps : mpsList )
				result.append(mps.toString());	
		} catch (java.io.IOException e) {
			e.printStackTrace();
		}
		return result.toString();
	}
}

class FRG extends ToString {
	String mid;
	String nm;
	byte [] seq;
	int lclr;
	int rclr;
}

class MPS {
	char typ;
	String mid;
	int lpos;
	int rpos;
	int dln; //gapno
	int del[]; //gaps
	public String toString() {
		StringBuilder result = new StringBuilder();
	     	try {
			FRG frg = Ca2ta.frgHash.get(mid);
			String seqDecompressed = Compression.decompress(frg.seq);
			System.out.println("Compressed size: " + frg.seq.length + ", Decompressed size: " + seqDecompressed.getBytes().length);
			result.append("#" + frg.nm + "(" + lpos + ")" + " [] " 
				+ seqDecompressed.length() 
				+ " bases, 00000000 checksum." 
				+ " {" + frg.lclr + "," + frg.rclr + "}" 
				+ " <" + lpos + "," + rpos +">" + "\n");
			result.append(Ca2ta.formatSeq(seqDecompressed));
		} catch (java.io.IOException e) {
			e.printStackTrace();
		}
		return result.toString();
	}
}
public class Ca2ta {
   static final int SEQ_OUTPUT_SIZE = 60;
   static Pattern acc = Pattern.compile("^acc:\\((\\d+)");
   static Pattern bracket = Pattern.compile("^\\{(\\w+)$");
   static Pattern field = Pattern.compile("^(\\w\\w\\w):(.*)$");

   static Map<String,FRG> frgHash;
   static Map<String,Contig> ccoHash;
   static String prefix;

   public static void main( String [] argv ) {
	if ( argv.length != 1 ) return;
	else   prefix = argv[0];

	//Step 0: Create clr file
	createCLR();
	
	//Step 1: Read prefix.clv file
	System.out.println("\nreadCLR");
	readClr();

	//Step 2: Read the FRG file
	System.out.println("\nreadFRG");
	readFrg();


	//Step 2: Read the noAFG file
	System.out.println("\nreadASM");
	readASM();

	//Step 3: Print the contig file
	System.out.println("\nprintContig");
	printContig();

   }
   
   static void createCLR() {
   	//First check if it already exists
		String clrFile = prefix.concat(".clr");
	
	if ( new File(clrFile).exists() )
		return;
   	
	try
        {   String cmd = "/usr/local/devel/ATG/moweis/CA_Latest/src/AS_RUN/arun/tools/asmToCLR.sh  " + prefix.concat(".asm");
	    System.out.println("Executing: " + cmd); //System.exit(1);
            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(cmd);
	    BufferedReader ls_in = new BufferedReader(
                                          new InputStreamReader(proc.getInputStream()));
	    FileWriter output = new FileWriter(clrFile);
	    BufferedWriter bufWrite = new BufferedWriter(output);
	    String ls_str;

	    try {
		while ((ls_str = ls_in.readLine()) != null) {
		    output.write(ls_str+"\n");
		}
	    } catch (IOException e) {
		System.exit(0);
	    }	    
	    ls_in.close();
	    bufWrite.close();
        } catch (Throwable t)
        {
            t.printStackTrace();
        }   
   }
   
   static void readClr() {
	String line = null;    // String that holds current file line
	try {
		FileReader input = new FileReader(prefix.concat(".clr"));

		BufferedReader bufRead = new BufferedReader(input);
		int count = 0;  // Line number of count 

		line = bufRead.readLine();

		frgHash = new HashMap<String,FRG>();

		while (line != null){
		  String[] values = line.split("\\s");
		  FRG frg = new FRG();
		  frg.mid = values[0];
		  frg.lclr = Integer.parseInt(values[1]);
		  frg.rclr = Integer.parseInt(values[2]);
		  frgHash.put(frg.mid,frg);
		  line = bufRead.readLine();
		  count++;
		  if ( count % 100000 == 0 ) System.out.println("CLR: " + count);

		}
		System.out.println("CLR Final Read: " + count);
		System.out.println("CLR Hash size: " + frgHash.size());
		bufRead.close();              
	}catch (ArrayIndexOutOfBoundsException e){
		System.out.println("Error reading line: " + line + "	\n");
		e.printStackTrace();
	}catch (IOException e){
		e.printStackTrace();
	}
   }

   static void readFrg() {
	try {
		FileReader input = new FileReader(prefix.concat(".frg"));

		BufferedReader bufRead = new BufferedReader(input);
   		String line = bufRead.readLine();      

        	int count = 0;
		
		while ( line != null ) {
			Matcher m = bracket.matcher(line);
			String recName = "";
			if ( m.lookingAt() ) {
				recName = m.group(1);
			}
			int brackets = 0;
			if (recName.equals("FRG") ) { 
				FRG frg = processFRG(bufRead);
				if ( frg != null )  {
					FRG getFRG = frgHash.get(frg.mid);
					if ( getFRG == null ) continue;
					getFRG.nm = frg.nm;
					getFRG.seq = frg.seq;
					count++;
			   	        if ( count % 100000 == 0 ) System.out.println("FRG: " + count);
				}
			} else {
				skipMsg(bufRead);
			}
			line = bufRead.readLine();
		}

		System.out.println("FRG Final Read: " + count);
		System.out.println("FRG Hash size: " + frgHash.size());
		bufRead.close();              	
	}catch (ArrayIndexOutOfBoundsException e){
		System.out.println("Usage: java ReadFile filename\n");
		e.printStackTrace();
	}catch (IOException e){
		e.printStackTrace();
	}catch (Exception e){
		e.printStackTrace();
	}
   }
   
   static void readASM() {
	try {
		FileReader input = new FileReader(prefix.concat(".asm"));
		BufferedReader bufRead = new BufferedReader(input);
		ccoHash = new HashMap<String,Contig>();
		getCARecord(bufRead);	
		bufRead.close();              
		System.out.println("Final CCO hash size: " + ccoHash.size());
	}catch (ArrayIndexOutOfBoundsException e){
		System.out.println("Usage: java ReadFile filename\n");
		e.printStackTrace();
	}catch (IOException e){
		e.printStackTrace();
	}catch (Exception e){
		e.printStackTrace();
	}   
   }
   
   static void printContig() {
		try {
			FileWriter output = new FileWriter(prefix.concat(".contig.BETA2"));
			BufferedWriter bufWrite = new BufferedWriter(output);
			writeContigs(bufWrite);	
			bufWrite.close();              
		}catch (ArrayIndexOutOfBoundsException e){
			System.out.println("Usage: java ReadFile filename\n");
			e.printStackTrace();
		}catch (IOException e){
			e.printStackTrace();
		}catch (Exception e){
			e.printStackTrace();
		}
   }   
   
   static void getCARecord ( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();
	int count = 0;
	while ( line != null ) {
		Matcher m = bracket.matcher(line);
		String recName = "";
		if ( m.lookingAt() ) {
			recName = m.group(1);
		}

		int brackets = 0;
		if (recName.equals("CCO") ) { 
			processContig(bufRead);
			count++;
			if ( count % 10000 == 0 ) System.out.println("CCO: " + count);
		} else {
			skipMsg(bufRead);
		}
		line = bufRead.readLine();
	}        
   }
   
   static FRG processFRG ( BufferedReader bufRead ) throws Exception {
   
   	String line = bufRead.readLine();
	FRG frg = new FRG();

	while ( line != null && !line.equals("}")) {
		Matcher fieldMatch = field.matcher(line);
		String fieldName = "";
		String fieldValue;
		try {
		if ( fieldMatch.lookingAt() )  {
			fieldName = fieldMatch.group(1);
			fieldValue = fieldMatch.group(2);			
			if ( fieldName.equals("acc") ) {
				frg.mid = fieldValue;
			} else if ( fieldName.equals("src") ) {
				frg.nm = bufRead.readLine();
			} else if ( fieldName.equals("seq") ) {			
				frg.seq = Compression.compress(readSequence(bufRead));
			}
		}
		} catch (NumberFormatException e)  {
			System.out.println("NumberFormatException at line: " + line);
			e.printStackTrace();
		}
		line = bufRead.readLine();
	}
	return frg;   
   }
   
   static void processContig( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();

	Matcher accMatch = acc.matcher(line);
	accMatch.lookingAt();
	Contig cco = new Contig();
	cco.acc = accMatch.group(1);	
	cco.mpsList = new ArrayList<MPS>();
	while ( line != null && !line.equals("}")) {
		Matcher bracketMatch = bracket.matcher(line);
		Matcher fieldMatch = field.matcher(line);
		String recName = "";
		String fieldName = "";
		String fieldValue;
		if ( bracketMatch.lookingAt() )
			recName = bracketMatch.group(1);
			if ( recName.equals("VAR") || recName.equals("UPS")) {
				skipMsg(bufRead);
			}
			else if ( recName.equals("MPS") ) {
				MPS returnedMPS = readMPS(bufRead);
				if ( returnedMPS != null ) {
					cco.mpsList.add(returnedMPS);
				}
			}
		else if ( fieldMatch.lookingAt() )  {
			fieldName = fieldMatch.group(1);
			fieldValue = fieldMatch.group(2);			
			if ( fieldName.equals("len") ) {
				cco.len = Integer.parseInt(fieldValue);
			} else if ( fieldName.equals("cns") ) {
				cco.cns = Compression.compress(readSequence(bufRead));
			} else if ( fieldName.equals("npc") ) {			
				cco.npc = Integer.parseInt(fieldValue);
			}
		} 
		line = bufRead.readLine();
	}
	ccoHash.put(cco.acc,cco);   
   }
   
   static MPS readMPS ( BufferedReader bufRead) throws Exception {
   	String line = bufRead.readLine();
	
   	MPS mps = new MPS();
	while ( line != null && !line.equals("}")) {
		Matcher fieldMatch = field.matcher(line);
		String fieldName = "";
		String fieldValue;
		if ( fieldMatch.lookingAt() )  {
			fieldName = fieldMatch.group(1);
			fieldValue = fieldMatch.group(2);			
			if ( fieldName.equals("typ") ) {
				if ( !fieldValue.equals("R") ) {
					return null;
				}
				mps.typ = fieldValue.charAt(0);
			} else if ( fieldName.equals("mid") ) {
				mps.mid = fieldValue;
			} else if ( fieldName.equals("pos") ) {			
				String [] posArray = fieldValue.split(",");
				mps.lpos = Integer.parseInt(posArray[0]);
				mps.rpos = Integer.parseInt(posArray[1]);
			} else if ( fieldName.equals("dln") ) {
				mps.dln = Integer.parseInt(fieldValue);
			} else if ( fieldName.equals("del") ) {
				line = bufRead.readLine();
				if ( line.equals("}") )  {
					return null;
				}
				String delStr = new String();
				while ( line != null && !line.equals("}") ) {
					delStr = delStr.concat(line).concat(" ");
					line = bufRead.readLine();
				}
				String [] delStrArray = delStr.split(" ");
				mps.del = new int[delStrArray.length];
				for( int i= 0; i < delStrArray.length ; i++ )
					mps.del[i] = Integer.parseInt(delStrArray[i]);
				if ( line.equals("}") )
					break;
			}
		} 
		line = bufRead.readLine();
	}
	return mps;
   }
   
   static void skipMsg ( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();
	int bracketCount = 1;
	while ( line != null) {
		if ( line.equals("}") ) {
			bracketCount--;
		}
		else {
			Matcher bracketMatch = bracket.matcher(line);
			if ( bracketMatch.lookingAt() )
				bracketCount++;		
		}
		
		if ( bracketCount == 0)
			return;
		line = bufRead.readLine();
	}   
   }
   
   
   static String readSequence ( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();
	String sequence = "";
	
	while ( line != null && !line.equals(".")) {
		sequence = sequence.concat(line);
		line = bufRead.readLine();
	}	
	
   	return sequence;
   }  
   
   static void writeContigs ( BufferedWriter bufWrite ) throws Exception {
	for (String contigId :  new TreeSet<String>(ccoHash.keySet()))
		bufWrite.write(ccoHash.get(contigId).toString());
   }
   
   static String formatSeq (String seq) {
	StringBuilder result = new StringBuilder();
		
   	if ( seq != null )
   		for (int i = 0 ; i < seq.length() ; i+= SEQ_OUTPUT_SIZE ) {
			int diff = seq.length() - i;
			if ( diff <= 0 ) break;
			int endpoint = ( diff < SEQ_OUTPUT_SIZE ) ? diff : SEQ_OUTPUT_SIZE;
			result.append(seq.substring(i,i+endpoint) + "\n");
		}
	return result.toString();
   }
}