# Script to get list of users in and out of the big summer teams

# split to teams and everyone else
teams = readRDS("data/BC_Parks_Big_Teams/BCParks_obs_unique_with_BigTeams_flagged.RDS") 
unique(teams$observer_type)

users_logins = teams |>
  group_by(observer_type, user_login) |>
  summarise("n" = n())

team = users_logins$user_login[which(users_logins$observer_type == "BC_BigTeam")]
visitors = users_logins$user_login[which(users_logins$observer_type != "BC_BigTeam")]
write.csv(data.frame("user_login" = team), "outputs/bigteam_userlogins.csv")