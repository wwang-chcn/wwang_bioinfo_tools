Storage()
{
 	user=$1
 	useradd -s /bin/bash -d /mnt/Storage/home/$user -m $user -g zhanglab
 	echo $user:123 > temp.passwd
 	/usr/sbin/chpasswd < temp.passwd
 	chage -d 0 $user
 
 	chown -R $user:zhanglab /mnt/Storage/home/$user
 	chmod 775 /mnt/Storage/home/$user

	usermod -a -G docker $user
 	
 	echo "$user is created"
 }

removeuser()
{
	user=$1
	userdel $user
	echo "$user is deleted"
	while true; do
		read -p "Do you wish to delete ${user}'s HOME as well? [Yes|NO]" yn
		case $yn in
			Yes ) rm -rf /mnt/Storage/home/$user && echo "${user}'s HOME is deleted"; break;;
			No ) exit;;
			* ) echo "Please answer Yes or No.";;
		esac
	done
}

chgroups()
{
	user=$1
	usermod -a -G sudo ${user}
}

# ----- Jan-16-2019 -----
# addgroup XXX_lab
# addgroup docker

# Storage admin_user1
# Storage admin_user2

# chgroups admin_user1
# chgroups admin_user2

# Storage normal_user1
# Storage normal_user2
# Storage normal_user3
# Storage normal_user4
# removeuser normal_user4
