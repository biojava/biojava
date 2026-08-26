package org.biojava.nbio.core.util;

import static org.biojava.nbio.core.util.FileDownloadUtils.getFileExtension;
import static org.biojava.nbio.core.util.FileDownloadUtils.getFilePrefix;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

class FileDownloadUtilsTest {

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

    		// files.wwpdb.org returns the content MD5 as the ETag, so the default
    		// ETag policy records a real checksum without a separate hash URL.
    		assertTrue(hashFile.exists(), "no hash file was derived from the ETag");
    		assertTrue(FileDownloadUtils.validateFile(destFile), "correctly downloaded file failed hash validation");

    		PrintStream temp2 = new PrintStream(hashFile);
    		temp2.print("ABCD"); // not a digest of any supported length
    		temp2.close();
    		// An unreadable sidecar must not condemn an otherwise good download.
    		assertTrue(FileDownloadUtils.validateFile(destFile),
    				"a malformed hash file should be ignored, not treated as a mismatch");

    		PrintStream temp3 = new PrintStream(hashFile);
    		temp3.print("00000000000000000000000000000000"); // well-formed but wrong MD5
    		temp3.close();
    		assertFalse(FileDownloadUtils.validateFile(destFile),
    				"file not detected to be invalid although hash value is wrong.");
    		System.out.println("Just ignore the previous warning. It is expected.");

    		destFile.delete();
    		sizeFile.delete();
    		hashFile.delete();
    	}
    }

    @Nested
    class HttpStatus {

        @Test
        void notFoundThrowsAndLeavesNothingBehind() throws IOException {
            // A path that is guaranteed absent from the wwPDB archive.
            URL missing = new URL("https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/zz/zzzz.cif.gz");
            File dest = new File(System.getProperty("java.io.tmpdir"), "bj-missing.cif.gz");
            File sizeFile = new File(dest.getParentFile(), dest.getName() + ".size");
            dest.delete();
            sizeFile.delete();

            HttpStatusException e = assertThrows(HttpStatusException.class,
                    () -> FileDownloadUtils.downloadFile(missing, dest));
            assertEquals(404, e.getStatusCode());
            assertTrue(e.isNotFound());
            assertFalse(dest.exists(), "a 404 body must never be written to the destination");

            // ... and no validation metadata may be recorded for it either, or the
            // cached error page would later pass validation.
            FileDownloadUtils.createValidationFiles(missing, dest, null, FileDownloadUtils.Hash.UNKNOWN);
            assertFalse(sizeFile.exists(), "no size file should be written for a 404 response");
        }
    }

    @Nested
    class Hashing {

        private File writeTemp(String name, byte[] content) throws IOException {
            File f = new File(System.getProperty("java.io.tmpdir"), name);
            Files.write(f.toPath(), content);
            f.deleteOnExit();
            return f;
        }

        @Test
        void digestsOfEmptyFileMatchKnownValues() throws IOException {
            File empty = writeTemp("bj-empty.bin", new byte[0]);
            assertEquals("d41d8cd98f00b204e9800998ecf8427e",
                    FileDownloadUtils.computeHash(empty, FileDownloadUtils.Hash.MD5));
            assertEquals("da39a3ee5e6b4b0d3255bfef95601890afd80709",
                    FileDownloadUtils.computeHash(empty, FileDownloadUtils.Hash.SHA1));
            assertEquals("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
                    FileDownloadUtils.computeHash(empty, FileDownloadUtils.Hash.SHA256));
        }

        @Test
        void digestOfKnownContent() throws IOException {
            File abc = writeTemp("bj-abc.bin", "abc".getBytes(StandardCharsets.UTF_8));
            assertEquals("900150983cd24fb0d6963f7d28e17f72",
                    FileDownloadUtils.computeHash(abc, FileDownloadUtils.Hash.MD5));
            assertTrue(FileDownloadUtils.verifyHash(abc, FileDownloadUtils.Hash.MD5,
                    "900150983CD24FB0D6963F7D28E17F72"), "comparison should be case-insensitive");
            assertFalse(FileDownloadUtils.verifyHash(abc, FileDownloadUtils.Hash.MD5,
                    "00000000000000000000000000000000"));
        }

        @Test
        void algorithmNames() {
            assertEquals("MD5", FileDownloadUtils.getAlgorithmName(FileDownloadUtils.Hash.MD5));
            assertEquals("SHA-1", FileDownloadUtils.getAlgorithmName(FileDownloadUtils.Hash.SHA1));
            assertEquals("SHA-256", FileDownloadUtils.getAlgorithmName(FileDownloadUtils.Hash.SHA256));
            assertThrows(IllegalArgumentException.class,
                    () -> FileDownloadUtils.getAlgorithmName(FileDownloadUtils.Hash.UNKNOWN));
        }
    }

    @Nested
    class HashFileParsing {

        private static final String MD5 = "900150983cd24fb0d6963f7d28e17f72";

        private String parse(String content) throws IOException {
            File f = new File(System.getProperty("java.io.tmpdir"), "bj-hashfile.txt");
            Files.write(f.toPath(), content.getBytes(StandardCharsets.UTF_8));
            f.deleteOnExit();
            return FileDownloadUtils.parseHashFile(f);
        }

        @Test
        void bareHex() throws IOException {
            assertEquals(MD5, parse(MD5));
            assertEquals(MD5, parse(MD5 + "\n"));
        }

        @Test
        void uppercaseHexIsKeptVerbatim() throws IOException {
            assertEquals(MD5.toUpperCase(), parse(MD5.toUpperCase()));
        }

        @Test
        void coreutilsLayouts() throws IOException {
            assertEquals(MD5, parse(MD5 + "  somefile.cif.gz\n"));
            assertEquals(MD5, parse(MD5 + " *somefile.cif.gz\n"));
        }

        @Test
        void bsdLayout() throws IOException {
            assertEquals(MD5, parse("MD5 (somefile.cif.gz) = " + MD5 + "\n"));
        }

        @Test
        void blankLeadingLinesAreSkipped() throws IOException {
            assertEquals(MD5, parse("\n   \n" + MD5 + "\n"));
        }

        @Test
        void garbageYieldsNull() throws IOException {
            assertNull(parse("not a hash at all\n"));
            assertNull(parse("ABCD\n"));
            assertNull(parse(""));
        }
    }

    @Nested
    class ETagParsing {

        @Test
        void wwpdbStyleMd5IsRecognised() {
            assertEquals(FileDownloadUtils.Hash.MD5,
                    FileDownloadUtils.hashFromETag("\"f99fb9d964e1e1c22f2ea559ac5745cf\""));
            assertEquals("f99fb9d964e1e1c22f2ea559ac5745cf",
                    FileDownloadUtils.normalizeETag("\"f99fb9d964e1e1c22f2ea559ac5745cf\""));
        }

        @Test
        void sha1AndSha256LengthsAreRecognised() {
            assertEquals(FileDownloadUtils.Hash.SHA1,
                    FileDownloadUtils.hashFromETag("da39a3ee5e6b4b0d3255bfef95601890afd80709"));
            assertEquals(FileDownloadUtils.Hash.SHA256,
                    FileDownloadUtils.hashFromETag(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"));
        }

        @Test
        void ebiStyleTimeSizeETagIsNotMistakenForADigest() {
            // nginx and Apache emit <hex-mtime>-<hex-size>; the dash keeps it out of
            // the hex-only pattern, which is what stops us recording a bogus checksum.
            assertEquals(FileDownloadUtils.Hash.UNKNOWN, FileDownloadUtils.hashFromETag("\"67910788-1009e0\""));
            assertEquals(FileDownloadUtils.Hash.UNKNOWN, FileDownloadUtils.hashFromETag("\"14c95dd-51eedb9922b40\""));
        }

        @Test
        void weakAndMissingETags() {
            assertEquals(FileDownloadUtils.Hash.UNKNOWN, FileDownloadUtils.hashFromETag("W/\"abc\""));
            assertEquals(FileDownloadUtils.Hash.UNKNOWN, FileDownloadUtils.hashFromETag(null));
            assertEquals(FileDownloadUtils.Hash.UNKNOWN, FileDownloadUtils.hashFromETag("  "));
            assertNull(FileDownloadUtils.normalizeETag(null));
            assertNull(FileDownloadUtils.normalizeETag("\"\""));
        }
    }

    @Nested
    class ValidateFile {

        @Test
        void bareRelativeNameDoesNotThrow() {
            // getParentFile() is null here; this used to be a NullPointerException.
            assertTrue(FileDownloadUtils.validateFile(new File("no-such-file-in-cwd.cif")));
        }

        @Test
        void emptySizeFileIsIgnoredRatherThanThrowing() throws IOException {
            File dir = Files.createTempDirectory("bj-validate").toFile();
            try {
                File data = new File(dir, "data.bin");
                Files.write(data.toPath(), "hello".getBytes(StandardCharsets.UTF_8));
                Files.write(new File(dir, "data.bin.size").toPath(), new byte[0]);
                assertTrue(FileDownloadUtils.validateFile(data));
            } finally {
                FileDownloadUtils.deleteDirectory(dir.toPath());
            }
        }

        @Test
        void sizeMismatchIsDetected() throws IOException {
            File dir = Files.createTempDirectory("bj-validate").toFile();
            try {
                File data = new File(dir, "data.bin");
                Files.write(data.toPath(), "hello".getBytes(StandardCharsets.UTF_8));
                Files.write(new File(dir, "data.bin.size").toPath(), "99".getBytes(StandardCharsets.UTF_8));
                assertFalse(FileDownloadUtils.validateFile(data));
            } finally {
                FileDownloadUtils.deleteDirectory(dir.toPath());
            }
        }

        @Test
        void everyHashSidecarIsChecked() throws IOException {
            File dir = Files.createTempDirectory("bj-validate").toFile();
            try {
                File data = new File(dir, "data.bin");
                Files.write(data.toPath(), "abc".getBytes(StandardCharsets.UTF_8));
                // correct MD5, wrong SHA1: the second sidecar must still be caught
                Files.write(new File(dir, "data.bin.hash_MD5").toPath(),
                        "900150983cd24fb0d6963f7d28e17f72".getBytes(StandardCharsets.UTF_8));
                Files.write(new File(dir, "data.bin.hash_SHA1").toPath(),
                        "0000000000000000000000000000000000000000".getBytes(StandardCharsets.UTF_8));
                assertFalse(FileDownloadUtils.validateFile(data));
            } finally {
                FileDownloadUtils.deleteDirectory(dir.toPath());
            }
        }

        @Test
        void writeHashFileRoundTrip() throws IOException {
            File dir = Files.createTempDirectory("bj-validate").toFile();
            try {
                File data = new File(dir, "data.bin");
                Files.write(data.toPath(), "abc".getBytes(StandardCharsets.UTF_8));
                FileDownloadUtils.writeHashFile(data, FileDownloadUtils.Hash.MD5,
                        FileDownloadUtils.computeHash(data, FileDownloadUtils.Hash.MD5));
                assertTrue(new File(dir, "data.bin.hash_MD5").exists());
                assertTrue(FileDownloadUtils.validateFile(data));
            } finally {
                FileDownloadUtils.deleteDirectory(dir.toPath());
            }
        }
    }
}
