/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of either the BSD licence or the GNU Lesser General
 * Public Licence.  These should be distributed with the code. 
 * If you do not have copies see:
 *
 *      http://www.opensource.org/licenses/bsd-license.php
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
 
package org.biojava.utils;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.StringReader;
import java.io.Writer;
import java.util.StringTokenizer;


/**
 * Convenience methods for running external processes.  This class
 * offers wrappers around the <code>java.lang.Process</code> API,
 * but hides away the details of managing threads and process I/O.
 *
 * <h3>Example</h3>
 *
 * <pre>
 * StringWriter out = new StringWriter();
 * ProcessTools.exec(
 *          new String[] {"/usr/bin/wc", "-w"},
 *          "The quick brown fox jumps over the lazy dog",
 *          out,
 *          null
 * );
 * int numWords = Integer.parseInt(out.toString().trim());
 * </pre>
 *
 * @author Thomas Down
 * @author Francois Pepin
 * @since 1.4
 * @deprecated preferable to use org.biojava.utils.ExecRunner 
 * or the org.biojava.utils.process package.
 */

public class ProcessTools {

        /** Win NT/2K/MEPro require cmd.exe to run programs **/
      private static final String WINDOWS_NT_2000_COMMAND_1 = "cmd.exe";

      /** Win NT/2K/MEPro require the /C to specify what to run **/
      private static final String WINDOWS_NT_2000_COMMAND_2 = "/C";

      /** Win 9X/MEHome require cmd.exe to run programs **/
      private static final String WINDOWS_9X_ME_COMMAND_1 = "command.exe";

      /** Win 9X/MEHome require the /C to specify what to run **/
      private static final String WINDOWS_9X_ME_COMMAND_2 = "/C";

    // Dummy constructor since this is just a tools class  
      
    private ProcessTools() {
    }
      
   /**
     * Execute the specified command and wait for it to return.
     *
     * @param args the command line to execute.
     * @param input data to present to the process' standard input, or <code>null</code> if the process does not require input.
     * @param stdout a <code>Writer</code> which will be filled with data from the process' output stream, or <code>null</code> to ignore output.
     * @param stderr a <code>Writer</code> which will be filled with data from the process' error stream, or <code>null</code> to ignore output.
     * @return the process' return code.
     */
     
    public static int exec(
        String[] args,
        Reader input,
        Writer stdout,
        Writer stderr
    )
        throws IOException
    {
      try {
          return exec(args, null, null, input, stdout, stderr, 0L);
        } catch (ProcessTimeoutException ex) {
            throw new Error("Assertion failed: unexpected process timeout");
        }
    }
  
   /**
     * Execute the specified command and wait for it to return, or kill
     * it if the specified timeout expires first.
     *
     * @param args the command line to execute.
     * @param envp environment variables for the child process, or <code>null</code> to inherit the current set.
     * @param dir working directory for the child process, or <code>null</code> to inherit the current directory.
     * @param input data to present to the process' standard input, or <code>null</code> if the process does not require input.
     * @param stdout a <code>Writer</code> which will be filled with data from the process' output stream, or <code>null</code> to ignore output.
     * @param stderr a <code>Writer</code> which will be filled with data from the process' error stream, or <code>null</code> to ignore output.
     * @param timeout maximum run-time (in milliseconds) for the child process.  A value of 0 indicates no limit.
     * @return the process' return code.
     * @throws IOException if an error occurs while starting or communicating with the process
     * @throws ProcessTimeoutException if the child process was killed because its timeout had expired.
     */
     
  public static int exec(
        String[] args,
        String[] envp,
        File dir,
        Reader input,
        Writer stdout,
        Writer stderr,
        long timeout
    )
        throws IOException, ProcessTimeoutException
    {
        Process proc = Runtime.getRuntime().exec(args, envp, dir);
        CharPump outPump;
        
        CharPump inPump, errPump;
        
        if (input == null) {
            input = new StringReader("");
        }
        outPump = new PumpReaderToWriter(input, new OutputStreamWriter(proc.getOutputStream()));
        if (stdout == null) {
            inPump = new PumpReaderToNull(new InputStreamReader(proc.getInputStream()));
        } else {
            inPump = new PumpReaderToWriter(new InputStreamReader(proc.getInputStream()), stdout);
        }
        if (stderr == null) {
            errPump = new PumpReaderToNull(new InputStreamReader(proc.getErrorStream()));
        } else {
          errPump = new PumpReaderToWriter(new InputStreamReader(proc.getErrorStream()), stderr);
        }
        
        TimeBomb tb = null;
        if (timeout > 0) {
            tb = new TimeBomb(proc, timeout);
            tb.start();
        }
        
        outPump.start();
        inPump.start();
        errPump.start();
        
        int rc;
        try {
            rc = proc.waitFor();
            
            if (tb != null) {
                tb.defuse();
            }
            
            outPump.join();
            inPump.join();
            errPump.join();
        } catch (InterruptedException iex) {
            throw new IOException("Error waiting for process to complete");
        }
        
        if (tb != null && tb.fired()) {
            // ProcessTimeoutException trumps any IOExceptions, since odd exceptions
            // may be generated after the process is destroyed.
            throw new ProcessTimeoutException(rc);
        } else {
            checkException(outPump, "Output to child");
            checkException(inPump, "Input from child");
            checkException(errPump, "Errors from child");
            return rc;
        }
    }

  
  /**
   * Execute the specified command and wait for it to return. This is the
   * simplified version that tries to be nice and make your life easier. If
   * you know exactly what you want, you might want to use exec(String[],...)
   * instead.  
   *
   * @param command the command line to execute.
   * @param input data to present to the process' standard input, or
   * <code>null</code> if the process does not require input. 
   * @param stdout a <code>Writer</code> which will be filled with data
   * from the process' output stream, or <code>null</code> to ignore output. 
   * @param stderr a <code>Writer</code> which will be filled with data
   * from the process' error stream, or <code>null</code> to ignore output. 
   * @return the process' return code.
   * @throws IOException if an error occurs while starting or communicating
   * with the process 
   */
  public static int exec(
                         String command,
                         Reader input,
                         Writer stdout,
                         Writer stderr)
        throws IOException
  {
      try {
            return exec(command, null, null, input, stdout, stderr, 0L);
      } catch (ProcessTimeoutException ex) {
            throw new Error("Assertion failed: unexpected process timeout");
      }
  }
  
  /**
   * Execute the specified command and wait for it to return. This is the
   * simplified version that tries to be nice and make your life easier. If
   * you know exactly what you want, you might want to use exec(String[],...)
   * instead.  
   *
   * @param command the command line to execute.
   * @param input data to present to the process' standard input, or
   * <code>null</code> if the process does not require input. 
   * @param stdout a <code>Writer</code> which will be filled with data
   * from the process' output stream, or <code>null</code> to ignore output. 
   * @param stderr a <code>Writer</code> which will be filled with data
   * from the process' error stream, or <code>null</code> to ignore output. 
   * @param timeout maximum run-time (in milliseconds) for the child process.  A value of 0 indicates no limit.
   * @return the process' return code.
   * @throws IOException if an error occurs while starting or communicating
   * with the process 
   * @throws ProcessTimeoutException if the child process was killed because its timeout had expired.
   */
  public static int exec(
                         String command,
                         String[] envp,
                         File dir,
                         Reader input,
                         Writer stdout,
                         Writer stderr,
                         long timeout)
        throws IOException, ProcessTimeoutException
  {
    String[] cmd = null;
    // First determine the OS to build the right command string
    String osName = System.getProperty("os.name");
    	 if (osName.equals("Windows NT") || osName.equals("Windows 2000") ||
             osName.equals("Windows XP")) {
	     cmd = new String[3];
	     cmd[0] = WINDOWS_NT_2000_COMMAND_1;
	     cmd[1] = WINDOWS_NT_2000_COMMAND_2;
	     cmd[2] = command;
	 }
	 else if (
	     osName.equals("Windows 95")
		 || osName.equals("Windows 98")
		 || osName.equalsIgnoreCase("Windows ME")) {
	     cmd = new String[3];
	     cmd[0] = WINDOWS_9X_ME_COMMAND_1;
	     cmd[1] = WINDOWS_9X_ME_COMMAND_2;
	     cmd[2] = command;
	 }
	 else {
	     // Linux (and probably other *nixes) prefers to be called
	     // with each argument supplied separately, so we first
	     // Tokenize it across spaces as the boundary.
	     StringTokenizer st = new StringTokenizer(command, " ");
	     cmd = new String[st.countTokens()];
	     int token = 0;
	     while (st.hasMoreTokens()) {
		 String tokenString = st.nextToken();
		 //System.out.println(tokenString);
		 cmd[token++] = tokenString;
	     }
	 }
         return exec(cmd, envp, dir, input,stdout,stderr, timeout);
  }
  
    /**
     * Check the status of a Pump and re-throw any exception which may
     * have occured during its lifecycle
     *
  
    private static void checkException(Pump p, String msg)
        throws IOException
    {
        IOException ioe = p.getException();
        if (ioe != null) {
            throw new IOException("Exception processing " + msg);
        }
    }*/

  private static void checkException(CharPump p, String msg)
        throws IOException
    {
        IOException ioe = p.getException();
        if (ioe != null) {
            throw new IOException("Exception processing " + msg);
        }
    }
  
    /**
     * Thread which will kill the specified process if it is not defused before
     * the timeout expires.
     */
    
    private static class TimeBomb extends Thread {
        private volatile boolean ticking = false;
        private volatile boolean fired = false;
        private final long time; 
        private final Process victim;
        
        public TimeBomb(Process victim, long time) {
            this.time = time;
            this.victim = victim;
        }
        
        public void run() {
            synchronized(this) {
                try {
                    ticking = true;
                    wait(time);
                } catch (InterruptedException ex) {
                    System.err.println("Timebomb thread was interrupted -- this shouldn't happen");
                }
            }
            if (ticking) {
                // System.err.println("TimeBomb activated -- killing child process");
                fired = true;
                victim.destroy();
            }
        }
        
        public synchronized void defuse() {
            ticking = false;
            notifyAll();
        }
        
        public boolean fired() {
            return fired;
        }
    }
    
    /**
     * Base class for threads which pump bytes from a source to a sink.  Subclasses
     * must implement sourceData and sinkData.  They may also wish to override
     * shutdownHook, a dummy method which is called when the pump finishes (usually
     * because the end of the data source has been reached).
     *
    
    private static abstract class Pump extends Thread {
        private IOException err = null;
        
        /**
         * Read bytes of data from some data source.
         *
        
        protected abstract int sourceData(byte[] buf) throws IOException;
        
        /**
         * Write bytes of data to some data sink.
         *
        
        protected abstract void sinkData(byte[] buf, int len) throws IOException;
        
        /**
         * Perform any required tidying operations when the Pump's job has finished.
         *
        
        protected void shutdownHook() throws IOException {};
        
        public void run() {
            try {
                byte[] buf = new byte[256];
                int cnt;
                do {
                    cnt = sourceData(buf);
                    if (cnt > 0) {
                        sinkData(buf, cnt);
                    }
                } while (cnt >= 0);
                shutdownHook();
            } catch (IOException e) {
                this.err = e;
            }
        }
        
        public IOException getException() {
            return err;
        }
    }*/

  private static abstract class CharPump extends Thread {
        private IOException err = null;
        
        /**
         * Read bytes of data from some data source.
         */
        
        protected abstract int sourceData(char[] buf) throws IOException;
        
        /**
         * Write bytes of data to some data sink.
         */
        
        protected abstract void sinkData(char[] buf, int len) throws IOException;
        
        /**
         * Perform any required tidying operations when the Pump's job has finished.
         */
        
        protected void shutdownHook() throws IOException {};
        
        public void run() {
            try {
                char[] buf = new char[256];
                int cnt;
                do {
                    cnt = sourceData(buf);
                    if (cnt > 0) {
                        sinkData(buf, cnt);
                    }
                } while (cnt >= 0);
                shutdownHook();
            } catch (IOException e) {
                this.err = e;
            }
        }
        
        public IOException getException() {
            return err;
        }
    }


  /*  
    private static final class PumpStreamToStringBuffer extends Pump {
        private final InputStream is;
        private final StringBuffer sb;
        
        public PumpStreamToStringBuffer(InputStream is, StringBuffer sb) {
            super();
            this.is = is;
            this.sb = sb;
        }
        
        protected int sourceData(byte[] buf)
            throws IOException
        {
            return is.read(buf);
        }
        
        protected void sinkData(byte[] buf, int len) {
            sb.append(new String(buf, 0, len));
        }       
    }
    
    private static final class PumpStreamToNull extends Pump {
        private final InputStream is;
        
        public PumpStreamToNull(InputStream is) {
            super();
            this.is = is;
        }
        
        protected int sourceData(byte[] buf)
            throws IOException
        {
            return is.read(buf);
        }
        
        protected void sinkData(byte[] buf, int len) {
        }       
    }
  */
  private static final class PumpReaderToNull extends CharPump {
        private final Reader is;
        
    public PumpReaderToNull(Reader is) {
            super();
            this.is = is;
        }
        
        protected int sourceData(char[] buf)
            throws IOException
        {
            return is.read(buf);
        }
        
        protected void sinkData(char[] buf, int len) {
        }       
    }
  

  private static final class PumpReaderToWriter extends CharPump {
    private final Reader reader;
    private final Writer writer;

    public PumpReaderToWriter(Reader reader, Writer writer){
      this.reader=reader;
      this.writer=writer;
    }
    
    protected int sourceData(char[] buf)
      throws IOException
    {
      return reader.read(buf, 0, buf.length);
    }

    protected void sinkData(char[] buf, int len)
      throws IOException
    {
      writer.write(buf, 0, len);
      writer.flush();
    }

    protected void shutdownHook() 
      throws IOException
    {
      writer.close();
    }
    
  }
  
  /*
    private static final class PumpStreamToStream extends Pump {
        private final InputStream is;
        private final OutputStream os;
        
        public PumpStreamToStream(InputStream is, OutputStream os) {
            this.is = is;
            this.os = os;
        }
        
        protected int sourceData(byte[] buf)
            throws IOException
        {
            return is.read(buf);
        }
        
        protected void sinkData(byte[] buf, int len)
            throws IOException
        {
            os.write(buf, 0, len);
            os.flush();
        }
        
        protected void shutdownHook() 
            throws IOException
        {
            os.close();
        }
        }*/
}
