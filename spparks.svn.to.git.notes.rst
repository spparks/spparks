
===============
Cloning Spparks
===============

1st attempt I got this message
------------------------------
git svn clone svn+ssh://dev.sandia.gov/usr/local/svn/spparks -s  spparks.git.svn

	WARNING: --prefix is not given, defaulting to empty prefix.  This is probably
	not what you want! In order to stay compatible with regular remote-tracking
	refs, provide a prefix like --prefix=origin/ (remember the trailing slash),
	which will cause the SVN-tracking refs to be placed at refs/remotes/origin/*.
	NOTE: In Git v2.0, the default prefix will change from empty to 'origin/'.
	Initialized empty Git repository in
	/home/jamitch/jaks.git/spparks.git.svn/.git/


1st attempt I got this message
------------------------------

git svn clone svn+ssh://dev.sandia.gov/usr/local/svn/spparks -s --prefix=spparks.svn/ spparks.git.svn

.... bunch of svn commits scroll by; then the last is r940; then the message below.

	Auto packing the repository for optimum performance. You may also
	run "git gc" manually. See "git help gc" for more information.
	Counting objects: 7922, done.
	Delta compression using up to 16 threads.
	Compressing objects: 100% (7862/7862), done.
	Writing objects: 100% (7922/7922), done.
	Total 7922 (delta 6624), reused 0 (delta 0)
	Checked out HEAD:
	  svn+ssh://dev.sandia.gov/usr/local/svn/spparks/trunk r940
