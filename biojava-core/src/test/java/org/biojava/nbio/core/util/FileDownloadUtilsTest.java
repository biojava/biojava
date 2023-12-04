package org.biojava.nbio.core.util;

import static org.biojava.nbio.core.util.FileDownloadUtils.getFileExtension;
import static org.biojava.nbio.core.util.FileDownloadUtils.getFilePrefix;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URL;
import java.nio.file.Files;

import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

class FileDownloadUtilsTest {

    @Nested
    class FileCopy {   
        
        private File createSrcFile () throws IOException {
            byte [] toSave = new byte []{1,2,3,4,5};
            File src = Files.createTempFile("test", ".dat").toFile();
            try (FileOutputStream fos = new FileOutputStream(src);){
                fos.write(toSave);
            }
            return src;
        }

        @Test
        void copyFile() throws IOException {
            File src = createSrcFile();
            //sanity check
            assertEquals(5, src.length());
            File dest = Files.createTempFile("dest", ".dat").toFile();
            assertEquals(0, dest.length());
            FileDownloadUtils.copy(src, dest);
            assertEquals(5, dest.length());

            //original is unaffected
            assertEquals(5, src.length());
            
            // bytes are identical
            try (FileInputStream fis1 = new FileInputStream(src);
                FileInputStream fis2 =  new FileInputStream(dest)) {
                int b = -1;
                while (( b = fis1.read()) != -1) {
                    int b2 = fis2.read();
                    assertEquals (b, b2);
                }
            }
        }
    }

    @Nested
    class FileExtension {
        @Test
        void getExtensionHappyCase(){
            File someFile = new File("sequence.fasta");
            assertEquals("fasta", getFileExtension(someFile));
        }

        @Test
        void lastSuffixOnlyReturned(){
            File someFile = new File("sequence.1.a.fasta");
            assertEquals("fasta", getFileExtension(someFile));
        }
        
        @Test
        void fileNameEndingInDotReturnsEmptyString(){
            File someFile = new File("noExtension.");
            assertEquals("", getFileExtension(someFile));
        }

        @Test
        void hiddenFile(){
            File someFile = new File(".m2");
            assertEquals("m2", getFileExtension(someFile));
        }

        @Test
        void noExtension(){
            File someFile = new File("nameOnly");
            assertEquals("nameOnly", getFileExtension(someFile));
        }        
    }

    @Nested
    class GetFilePrefix{
        @Test
        void standardFileName(){
            File someFile = new File("sequence.fasta");
            assertEquals("sequence", getFilePrefix(someFile));
        }
        @Test
        void prefixIsUpToFirstDot(){
            File someFile = new File("sequence.1.2.fasta");
            assertEquals("sequence", getFilePrefix(someFile));
        }

        @Test
        void noExtension(){
            File someFile = new File("nameOnly");
            assertEquals("nameOnly", getFilePrefix(someFile));
        }

        @Test
        void hiddenFile(){
            File someFile = new File(".m2");
            assertEquals("", getFilePrefix(someFile));
        }
    }

    @Nested
    class ToUnixPath {
        @Test
        void windowsToUnixAddsTrailingSlash(){
            String winPath = "C:\\a\\b\\c";
            assertEquals("C:/a/b/c/", FileDownloadUtils.toUnixPath(winPath));
        }
        @Test
        void unixPathReturnedUnchanged(){
            String path = "/a/b/c/";
            assertEquals(path, FileDownloadUtils.toUnixPath(path));
        }
    }

    @Nested
    class ExpandUserHome {
        String currUserHome = System.getProperty("user.home");
        @Test
        void minimalPath (){
        	String path="~";
        	assertEquals(currUserHome, FileDownloadUtils.expandUserHome(path));
        }
        @Test
        void simplePath (){
            String path="~/sequence.gb";
            assertEquals(currUserHome+File.separator+"sequence.gb", FileDownloadUtils.expandUserHome(path));
        }
        @Test
        void nestedPath (){
            String path="~/a/b/c/sequence.gb";
            assertEquals(currUserHome+File.separator
            		+ "a" + File.separator 
            		+ "b" + File.separator 
            		+ "c" + File.separator 
            		+ "sequence.gb", 
            		FileDownloadUtils.expandUserHome(path));
        }  
    }

    @Nested
    class URLMethods {
        final String availableUrl = "https://www.google.com";

        @Test
        void pingGoogleOK(){
            assertTrue(FileDownloadUtils.ping(availableUrl, 1000));
        }

        @Test
        void pingNonExistentFalse(){
            assertFalse(FileDownloadUtils.ping("https://non-existent.biojava", 1));
        }
    }
    @Nested
    class DeleteDirectory {

        private File createDirectoryTree () throws IOException {

            File tmpdir = Files.createTempDirectory("tmpDirPrefix").toFile();
            File child1 = new File(tmpdir, "a");
            File child2 = new File(child1, "b");
            File child3 = new File(child2, "c");
            File f = new File(child3, "seq.fa");
            child3.mkdirs();
            f.createNewFile();
            return tmpdir;
        }

        @Test
        void deleteFolderTree() throws IOException{
            File toDelete = createDirectoryTree();
            assertTrue(toDelete.exists());

            FileDownloadUtils.deleteDirectory(toDelete.getAbsolutePath());
            assertFalse(toDelete.exists());
        }
    }
    
    @Nested
    class CreateValidationFiles{
    	
    	@Test
    	void testValidationFiles() throws IOException{
    		URL sourceUrl = new URL("https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/45/145d.cif.gz");
    		File destFile = new File(System.getProperty("java.io.tmpdir"), "145d.cif.gz");
    		File sizeFile = new File(destFile.getParentFile(), destFile.getName()+".size");
    		File hashFile = new File(destFile.getParentFile(), destFile.getName()+".hash_MD5");
    		System.out.println(destFile.getAbsolutePath());
    		destFile.delete();
    		sizeFile.delete();
    		hashFile.delete();
    		assertFalse(destFile.exists(), "couldn't delete dest file");
    		assertFalse(sizeFile.exists(), "couldn't delete size file");
    		assertFalse(hashFile.exists(), "couldn't delete hash file");
    		
    		FileDownloadUtils.downloadFile(sourceUrl, destFile);
    		assertTrue(destFile.exists(), "couldn't create dest file");

    		assertTrue(FileDownloadUtils.validateFile(destFile), "file detected to be invalid although there are no validation files");

    		PrintStream temp1 = new PrintStream(sizeFile);
    		temp1.print(15); // some wrong size value
    		temp1.close();
    		assertFalse(FileDownloadUtils.validateFile(destFile), "file not detected to be invalid although size value is wrong.");
    		System.out.println("Just ignore the previous warning. It is expected.");
    		
    		FileDownloadUtils.createValidationFiles(sourceUrl, destFile, null, FileDownloadUtils.Hash.UNKNOWN);
    		assertTrue(sizeFile.exists(), "couldn't create size file");
    		assertTrue(FileDownloadUtils.validateFile(destFile), "file not detected to be invalid although there is correct size validation file");

    		PrintStream temp2 = new PrintStream(hashFile);
    		temp2.print("ABCD"); // some wrong hash value
    		temp2.close();
    		//This is not yet implemented. I am using this test for documentation purpose.
    		assertThrows(UnsupportedOperationException.class, 
    				() -> FileDownloadUtils.validateFile(destFile), 
    				"file not detected to be invalid although hash value is wrong.");
    		
    		destFile.delete();
    		sizeFile.delete();
    		hashFile.delete();
    	}
    }
}
