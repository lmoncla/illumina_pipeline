## Bitbucket and git tutorial  
This guide is designed to teach you how to setup git and download the scripts to begin using the TCF Lab sequencing pipeline. This is just a primer, [bitbucket](www.bitbucket.org) and [github](www.github.com) have similar guides that are in much greater detail. Consult those if you have any issues with the following process.

#### Step 1: Verifying git installed
Open terminal and type `git --version` and hit enter.  

If you get `-bash: git: command not found` please go [here](https://www.atlassian.com/git/tutorials/install-git) to and follow the instructions for your operating system to install git.  

#### Step 2a: Getting the scripts
If you need to setup an account for connecting to and interacting with the TCF lab repository please goto step 2b.  

If you only need access to the sequencing pipeline do the following.

1. In terminal navigate to the directory where you want to scripts to be.

2. Enter the following command `git clone https://bitbucket.org/tcf-lab/illumina_pipeline.git`

#### Step 2b: Setup an account at [bitbucket](www.bitbucket.org)
First create an account at [bitbucket.org](www.bitbucket.org) and have a member of the lab add you to the TCF team by providing them with your username or email that you used to create you bitbucket account.

#### Step 3: Configuring git
If you have already used git on your machine skip this step. If your are unsure you can use the `git config --global --edit` to check your configuration.  
_Note: public repositories you create will have this information._

1. Use the `git config --global user.name <name>` and replace `<name>` with your name.

2. Use the `git config --global user.email <email>` and replace `<email>` with your email.

3. Use the `git config --system core.editor nano` to set nano as your default editor for comments.

#### Step 4: Setting up git to work with SSH
We need to set up your machine to work securely with bitbucket. We will do this using a network communication protocol known as SSH.

1. In terminal run the following command to generate a key pair `ssh-keygen`.  
_If you get an_ `-bash: ssh-keygen: command not found` _you will need to install SSH._  
The command will ask you for a passphrase. It is recommended to use a passphrase, however you will be required to use the passphrase anytime you want to upload or download from bitbucket. So I will leave it up to you if you want to use one, if not just leave it blank and hit enter.

2. Use the command `cat ~/.ssh/id_rsa.pub` and copy  to clipboard or if you are using osx the command `pbcopy < ~/.ssh/id_rsa.pub` will automatically copy the key to the clipboard.

3. Login to [bitbucket](www.bitbucket.org) and click on your avatar at the top right of the screen. Find "Bitbucket Settings", in this window find SSH Keys under security. Click add Key.

4. Label this key relevant to the computer you are using, the key will authenticate your computer to interact with the bitbucket repository. If you decide to work on a different computer you will have to do this process again to be able to upload/download from the repository.

5. Paste the key copied from 2 into the key field and hit add key.

#### Step 5: Downloading the repository.
Git is a valuable tool for version control and maintaining a repository of digital work for the lab. It's an important skill to learn how to use. In the interest of time this guide will only show you how to clone the repository from bitbucket.

1. Navigate to the TCF team page on [bitbucket](www.bitbucket.org) and click on projects. You will see a folder for each lab member (ask if yours hasn't been created yet). Click on the public folder.

2. Click on the illumina_pipeline repository. Click "clone" and copy the command and execute it in terminal where you want the directory to be.
