/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
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

/*
 *    ExternalProcess.java
 */
package org.biojava.utils.process;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Enumeration;
import java.util.Properties;

import org.biojava.utils.SimpleThreadPool;
import org.biojava.utils.ThreadPool;

/**
 * Utility class to execute an external process and to handle 
 * the <code>STDOUT</code>, <code>STDERR</code> and <code>STDIN</code> streams
 * in multiple threads managed by a thread pool.
 * <p>This class is intended for applications that call an external program many 
 * times, e.g. in a loop, and that need high performance throughput, i.e.
 * the program's input and output should not be written to disk. The Java
 * {@link java.lang.Runtime#exec} methods requires the application to read/write
 * the external program's input and output streams in multiple threads. 
 * Otherwise the calling application may block. However, instantiating multiple 
 * threads for each call is extensive. On Linux systems there is also the 
 * problem that each Java thread is represented by a single process and the
 * number of processes is limited on Linux. Because the Java garbage collector
 * does not free the {@link java.lang.Thread} objects properly, an application
 * might run out of threads (indicated by a {@link java.lang.OutOfMemoryError} 
 * exception) after multiple iterations. Therefore, the 
 * <code>ExternalProcess</code> class uses a 
 * {@linkplain org.biojava.utils.ThreadPool thread pool}.</p>
 * <p>The simplest way to use this class is by calling the static methods
 * {@link #execute(String)} and 
 * {@link #execute(String, String, StringWriter, StringWriter)}. However, these
 * methods are not thread safe and no configuration is possible. In the former
 * case the program's input, output and error output is redirected to <code>STDIN</code>,
 * <code>STDOUT</code> and <code>STDERR</code> of the calling program. In the 
 * latter case input is provided as string and output and error output is 
 * written to {@link java.io.StringWriter} objects. The environment, i.e.
 * the current working directory and the environment variables, are inherited
 * from the calling process. In both cases, a static thread pool of size 
 * {@link #THREAD_POOL_SIZE} is used. The command that should be executed is 
 * provided as a string argument.</p>
 * <p>In scenarios where the environment has to be changed, the program input
 * is generated just in time, or the program's output is parsed just in time,
 * the use of an explicit instance of the <code>ExternalProcess</code> class
 * is recommended. This instance could be initialized with a custom thread pool.
 * Otherwise a {@link org.biojava.utils.SimpleThreadPool} of size 3 is used.
 * The input and output is managed by multithreaded 
 * {@linkplain org.biojava.utils.process.InputHandler input handler} and 
 * {@linkplain org.biojava.utils.process.OutputHandler output handler} objects.
 * There are four predefined handlers that read the program's input from a
 * {@link java.io.Reader} object or a {@link java.io.InputStream} object and
 * write the program's output to a {@link java.io.Writer} object or a
 * {@link java.io.OutputStream} object. These classes are called:
 * {@link org.biojava.utils.process.ReaderInputHandler}, 
 * {@link org.biojava.utils.process.SimpleInputHandler},
 * {@link org.biojava.utils.process.WriterOutputHandler} and
 * {@link org.biojava.utils.process.SimpleOutputHandler}. If no handlers are
 * specified the input and output is redirected to the standards streams of 
 * the calling process.</p>
 * <p>Before one of the methods {@link #execute()} or 
 * {@link #execute(Properties)} is called, the {@linkplain #setCommands(String)
 * commands} property should be set. One may include placeholders of the form
 * <code>%PARAM%</code> within the commands. If a 
 * {@link java.util.Properties} object is passed to the 
 * {@link #execute(Properties)} method, the placeholders are replaced by the 
 * particular property value. Therefore, the <code>Properties</code> object
 * must contain a key named <code>PARAM</code> (case doesn't matter). The 
 * environment for calling the external program can be configured using the
 * properties {@linkplain #setWorkingDirectory(File) workingDirectory} and
 * {@linkplain #setEnvironmentProperties(String[]) environmentProperties}.</p>
 * <p>Finally, the {@linkplain #setSleepTime(int) sleepTime} property can be
 * increased, in case the output handlers are not able to catch the whole
 * program's output within the given time. The default value is 
 * {@link #SLEEP_TIME} [in milliseconds].</p> 
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 * @see java.lang.Process
 */
public final class ExternalProcess {
    
    /* STATIC FIELDS */
    
    /**
     * Size of the thread pool for the static execute methods.
     */
    public static final int THREAD_POOL_SIZE = 9;
    
    /**
     * Number of milliseconds the execute method should pauses after the 
     * external process has finished the execution.
     */
    public static final int SLEEP_TIME = 10;
    
    /**
     * Static thread pool for static execute method. 
     */
    private static final ThreadPool THREAD_POOL = 
        new SimpleThreadPool(THREAD_POOL_SIZE, true);
    
    
    /* MAIN METHOD */
    
    /**
     * Runs an external program from the command line. The external process
     * inherits the environment variables and the current working directory from
     * the parent process.
     * @param args the path or the name of the external program and its command
     * line arguments
     */
    public static void main(/*@non_null@*/ String[] args) {
        try {
            System.exit(execute(joinCommands(args)));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /* STATIC METHODS */
    
    /**
     * Resolves the given command line by replacing all placeholder of the 
     * format <code>%NAME%</code> with the values from the given properties
     * for the corresponding keys of the format <code>NAME</code>.
     * @param commands the given command line
     * @param variables the placeholders or <code>null</code> if no resolvement
     * should be performed
     * @return the new command line
     * @throws NullPointerException if <code>commands</code> is 
     * <code>null</code>.
     */
    public static String resolveCommands(/*@non_null@*/ String commands, 
            Properties variables) throws NullPointerException {
        if (commands == null) {
            throw new NullPointerException("commands is null.");
        }
        if (variables != null && !variables.isEmpty()) {
            Enumeration keys = variables.keys();
            while (keys.hasMoreElements()) {
                String key = (String) keys.nextElement();
                String value = variables.getProperty(key);
                // TODO quotting
                commands = commands.replaceAll("%" + key + "%", value);
            }
        }
        return commands;        
    }

    /**
     * Executes an external program. The working directory and the environment 
     * variables are inherited from the parent process. The program input is
     * read from <code>STDIN</code>, the program output is written to 
     * <code>STDOUT</code> and the program error output is written to
     * <code>STDERR</code>.
     * <p><b>Note:</b> This method is not thread-safe.</p>
     * @param commands the command line including the path
     * or the name of the external program and its command line arguments
     * @return the exit code from the external program
     * @throws SecurityException if a security manager exists and its 
     * <code>checkExec</code> method doesn't allow creation of a subprocess.
     * @throws IOException if an I/O error occurs.
     * @throws NullPointerException if <code>commands</code> is 
     * <code>null</code>.
     * @throws IllegalArgumentException if <code>commandList</code> is empty.
     * @throws InterruptedException if the current thread is 
     * {@link Thread#interrupt() interrupted} by another thread 
     * while it is waiting, then the wait is ended and an 
     * <code>InterruptedException</code> is thrown.
     */
    public static int execute(/*@non_null@*/ String commands) 
    throws IOException, InterruptedException, NullPointerException, 
    SecurityException, IllegalArgumentException {
        return execute(commands, null, null, null);
    }
    
    /**
     * Executes an external program. The working directory and the environment 
     * variables are inherited from the parent process.
     * <p><b>Note:</b> This method is not thread-safe.</p>
     * @param commands the command line including the path
     * or the name of the external program and its command line arguments
     * @param inputString the input for the external programm or
     * <code>null</code> if the input should be read from <code>STDIN</code>
     * @param outputString the output of the external programm or
     * <code>null</code> if the output should be written to <code>STDOUT</code>
     * @param errorString the error output of the external program or
     * <code>null</code> if the error output should be written to 
     * <code>STDERR</code>
     * @return the exit code from the external program
     * @throws SecurityException if a security manager exists and its 
     * <code>checkExec</code> method doesn't allow creation of a subprocess.
     * @throws IOException if an I/O error occurs.
     * @throws NullPointerException if <code>commandList</code> is 
     * <code>null</code>.
     * @throws IllegalArgumentException if <code>commandList</code> is empty.
     * @throws InterruptedException if the current thread is 
     * {@link Thread#interrupt() interrupted} by another thread 
     * while it is waiting, then the wait is ended and an 
     * <code>InterruptedException</code> is thrown.
     */
    public static int execute(/*@non_null@*/ String commands, 
            String inputString, StringWriter outputString, 
            StringWriter errorString) 
    throws IOException, InterruptedException, NullPointerException, 
    SecurityException, IllegalArgumentException {
        
        ExternalProcess ep = new ExternalProcess(THREAD_POOL);
      
        ep.setCommands(commands);
                  
        if (inputString != null) {
            ep.setInputHandler(new ReaderInputHandler(
                    new StringReader(inputString), "STDIN"));
        }
        if (outputString != null) {
            ep.setOutputHandler(new WriterOutputHandler(outputString, 
                    "STDOUT"));
        }
        if (errorString != null) {
            ep.setErrorHandler(new WriterOutputHandler(errorString, "STDERR"));
        }

	    return ep.execute();
    }    
    
    /**
     * Joins a command list to a single command string.
     * @param commandList the list of the command and its arguments
     * @return the joined command line
     * @throws NullPointerException if <code>commandList</code> is 
     * <code>null</code>.
     */
    public static String joinCommands(/*@non_null@*/ Object[] commandList) 
    throws NullPointerException {
        if (commandList == null) {
            throw new NullPointerException("commandList is null.");
        }
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < commandList.length; i++) {
            // TODO quoting
            if (i == commandList.length - 1) {
                sb.append(commandList[i].toString());
            } else {
                sb.append(commandList[i].toString() + " ");
            }
        }
        return sb.toString();
    }
    
    /* PRIVATE FIELDS */

	/**
	 * The command including the external program and its arguments.
	 */
	private String commands = null;
    
	/**
	 * The working directory of the external program.
	 */
	private File workingDirectory  = null;
	
    /**
     * The list of environment variables for the external program.
     */
    private String[] environmentProperties = null;
  
	/**
	 * The input handler for STDIN.
	 */
	private InputHandler inputHandler = null;
	
    /**
     * The output handler for STDOUT.
     */
    private OutputHandler outputHandler = null;
	
    /**
     * The output handler for STDERR.
     */
    private OutputHandler errorHandler = null;
	
	/**
	 * The thread pool for the input and output handlers.
	 */
	private ThreadPool threadPool = null;
    
    /**
     * Should threads be stopped when finalized?
     */
    private boolean stopThreads = true;
    
    /**
     * Number of seconds to wait for completion of stream handlers.
     */
    private int sleepTime = SLEEP_TIME;
        
    /* PUBLIC CONSTRUCTORS */
    
    /**
     * Initializes the external process. 
     */
    public ExternalProcess() {
        this(null);
    }

	/**
     * Initializes the external process. 
	 * @param threadPool a thread pool with at least three threads or 
     * <code>null</code> if the default thread pool should be used
	 */
	public ExternalProcess(ThreadPool threadPool) {
        if (threadPool == null) {
            this.threadPool = new SimpleThreadPool(3, true);
            stopThreads = true;
        } else {
            this.threadPool = threadPool;
            stopThreads = false;
        }
        setEnvironmentProperties(null);
        setWorkingDirectory(null);
        setInputHandler(null);
        setOutputHandler(null);
        setErrorHandler(null);
	}
    
    /* PUBLIC METHODS */
  
	/**
     * Executes the external process and waits for its termination.
	 * @return the exit code from the external process
     * @throws IllegalArgumentException if the command is empty
     * @throws SecurityException if a security manager exists and its 
     * <code>checkExec</code> method doesn't allow creation of a subprocess.
     * @throws IOException if an I/O error occurs.
     * @throws InterruptedException if the current thread is 
     * {@link Thread#interrupt() interrupted} by another thread 
     * while it is waiting, then the wait is ended and an 
     * <code>InterruptedException</code> is thrown.
	 */
	public /*@pure@*/ int execute() throws IOException, InterruptedException, 
    SecurityException, IllegalArgumentException {
	    return execute((Properties) null);
	}

    /**
     * Executes the external process and waits for its termination.
     * @param variables a list of key-value-pairs that should be used to replace
     * placeholders in the command line. May be <code>null</code>.
     * @return the exit code from the external process
     * @throws SecurityException if a security manager exists and its 
     * <code>checkExec</code> method doesn't allow creation of a subprocess.
     * @throws IllegalArgumentException if the command is empty
     * @throws IOException if an I/O error occurs.
     * @throws InterruptedException if the current thread is 
     * {@link Thread#interrupt() interrupted} by another thread 
     * while it is waiting, then the wait is ended and an 
     * <code>InterruptedException</code> is thrown.
     */
	public /*@pure@*/ int execute(Properties variables) 
    throws IOException, InterruptedException, SecurityException, 
    IllegalArgumentException {
    
		Runtime runtime = Runtime.getRuntime();
		String commands = resolveCommands(this.commands, variables);
		
		Process process = runtime.exec(commands, environmentProperties, 
                workingDirectory);
    	
	    OutputStream in = process.getOutputStream();
	    InputStream out = process.getInputStream();
	    InputStream err = process.getErrorStream();
		
		inputHandler.setOutput(in);
		outputHandler.setInput(out);
		errorHandler.setInput(err);
        
		threadPool.addRequest(inputHandler);
		threadPool.addRequest(outputHandler);
		threadPool.addRequest(errorHandler);
		
        // start input and output handlers
		Thread.yield();
	
		int exitCode = process.waitFor();
			
        // give stream handlers time to complete
        Thread.sleep(sleepTime);

        in.close();
		out.close();
		err.close();
		
		return exitCode;
	}

    /**
     * Gets the command line including the path or name of the external program
     * and its command line arguments.
     * @return the command line
     */ 
    public /*@pure non_null@*/ String getCommands() {
        return commands;
    }
    
    /**
     * Sets the command line including the path or name of the external program
     * and its command line arguments.
     * @param commands the command line
     * @throws NullPointerException if <code>commands</code> is 
     * <code>null</code>.
     */
    public void setCommands(/*@non_null@*/ String commands) 
    throws NullPointerException {
        if (commands == null) {
            throw new NullPointerException("commands is null.");
        }
        this.commands = commands;
    }
    
    /**
     * Gets environment variables for the external process.
     * @return a list of strings in the format
     * <code>name=value</code> or <code>null</code> if the environment variables
     * should be inherited from the parent process
     */
    public /*@pure@*/ String[] getEnvironmentProperties() {
        String[] environmentProperties = null;
        if (this.environmentProperties != null) {
            environmentProperties = 
                (String[]) this.environmentProperties.clone();
        }        
        return environmentProperties;
    }
    
    /**
     * Sets environment variables for the external process.
     * @param environmentProperties a list of strings in the format
     * <code>name=value</code> or <code>null</code> if the environment variables
     * should be inherited from the parent process
     */
    public void setEnvironmentProperties(String[] environmentProperties) {
        if (environmentProperties == null) {
            this.environmentProperties = null;
        } else {
            this.environmentProperties = 
                (String[]) environmentProperties.clone();
        }
    }

    /**
     * Gets the output error handler which is responsible for the standard error
     * output of the external process.
     * @return the error output handler
     */
    public /*@pure non_null@*/ OutputHandler getErrorHandler() {
        return errorHandler;
    }
    
    /**
     * Sets the output error handler which is responsible for the standard error
     * output of the external process.
     * @param errorHandler the error output handler or <code>null</code> if the
     * error output should be redirected to <code>STDERR</code>
     */
    public void setErrorHandler(OutputHandler errorHandler) {
	    if (errorHandler == null) {
	        this.errorHandler = new SimpleOutputHandler(System.err, "STDERR");
	    } else {
	        this.errorHandler = errorHandler;
	    }
    }

    /**
     * Gets the input handler which is responsible for the standard input
     * of the external process.
     * @return the input handler
     */
    public /*@pure non_null@*/ InputHandler getInputHandler() {
        return inputHandler;
    }
    
    /**
     * Sets the input handler which is responsible for the standard input
     * of the external process.
     * @param inputHandler the input handler or <code>null</code> if the
     * input should be read from <code>STDIN</code>
     */
    public void setInputHandler(InputHandler inputHandler) {
	    if (inputHandler == null) {
	        this.inputHandler = new SimpleInputHandler(System.in, "STDINS");
	    } else {
	        this.inputHandler = inputHandler;
	    }
    }
    
    /**
     * Gets the output handler which is responsible for the standard output
     * of the external process.
     * @return the output handler
     */
    public /*@pure non_null@*/ OutputHandler getOutputHandler() {
        return outputHandler;
    }
    
    /**
     * Sets the output handler which is responsible for the standard output
     * of the external process.
     * @param outputHandler the output handler or <code>null</code> if the
     * output should be redirected to <code>STDOUT</code>
     */
    public void setOutputHandler(OutputHandler outputHandler) {
	    if (outputHandler == null) {
	        this.outputHandler = new SimpleOutputHandler(System.out, "STDOUT");
	    } else {
	        this.outputHandler = outputHandler;
	    }
    }

    /**
     * Gets the thread pool which is used for the input and output handlers.
     * @return a thread pool with at least three threads
     */
    public /*@pure non_null@*/ ThreadPool threadPool() {
        return threadPool;
    }
    
    /**
     * Gets the working directory for the external process.
     * @return the working directory or <code>null</code> if it should be 
     * inherited from the parent process
     */
    public /*@pure@*/ File getWorkingDirectory() {
        return workingDirectory;
    }
    
    /**
     * Sets the working directory for the external process.
     * @param workingDirectory the working directory or <code>null</code> if it 
     * should be inherited from the parent process
     */
    public void setWorkingDirectory(File workingDirectory) {
        this.workingDirectory = workingDirectory;
    }
        
    /*@ ensures \result >= 0; @*/
    /**
     * Gets the number of milliseconds the {@linkplain #execute(Properties)}
     * method should pauses after the external process is terminated. This gives
     * the stream handlers the time to complete their work. 
     * @return time in milliseconds
     */
    public /*@pure@*/ int getSleepTime() {
        return sleepTime;
    }

    /*@ requires sleepTime >= 0; @*/
    /**
     * Sets the number of milliseconds the {@linkplain #execute(Properties)}
     * method should pauses after the external process is terminated. Increase
     * this value if the output handlers didn't catch the whole program output. 
     * @param sleepTime time in milliseconds
     * @throws IllegalArgumentException if <code>sleepTime</code> is negative. 
     */
    public void setSleepTime(int sleepTime) 
    throws IllegalArgumentException {
        if (sleepTime < 0.0) {
            throw new IllegalArgumentException(
                    "sleepTime must be zero or positive.");
        }
        this.sleepTime = sleepTime;
    }

    /* CLASS Object */

    /**
     * {@inheritDoc}
     */
    protected void finalize() {
        if (stopThreads) {
            threadPool.stopThreads();
        }
    }
}
